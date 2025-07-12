import yaml, os
from core.calculator import QCCalculator
from core.loader import MODEL, PredictorLoader
from core.utils import gen_atoms, gen_xyz, setup_logger, parse_args_into_config

def main(
    config: dict[str, dict]
    ):
    if config['general']['log_file']:
        log_file = config['general']['log_file']
    else:
        name, suffix = os.path.splitext(config['structure']['file_path'])
        log_file = name + '.log'
    logger = setup_logger('QC', log_file)

    try:
        logger.info('Loading calculator...')
        loader = PredictorLoader(model_size=config['model']['size'], device=config['model']['device'])
        calc = loader.calc
        logger.info('Done!')

        atoms = gen_atoms(
            config['structure']['file_path'],
            config['structure']['charge'],
            config['structure']['spin'])

        QCC = QCCalculator(atoms=atoms, calc=calc, logger=logger)
        logger.info(f'file path: {config["structure"]["file_path"]}')
        logger.info(f'Input orientation: \n{gen_xyz(atoms)}')
        logger.info(f'charge: {atoms.info["charge"]}')
        logger.info(f'spin: {atoms.info["spin"]}')
        logger.info(f'model: {MODEL[config["model"]["size"]]}')
        for task, args in config['task'].items():
            args: dict
            if args['run']:
                kwargs = {k: v for k, v in args.items() if k != 'run'}
                getattr(QCC, task)(**kwargs)

    except Exception as e:
        logger.error("Error occurred during execution:", exc_info=True)
        raise
    finally:
        if config['model'].get('device') == 'cuda':
            import torch
            torch.cuda.empty_cache()

if __name__ == '__main__':
    # TODO: 支持批量计算，类似 `qg09 -a`
    with open('config.yml', 'r', encoding='utf-8') as f:
        config: dict = yaml.full_load(f)
    config = parse_args_into_config(config)
    main(config)

