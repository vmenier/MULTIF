python runModel.py -l 1 -f general_gradients.cfg
mkdir Baseline
python train.py
rm -rf cv_data train.py train.dat shaowuML.py shaowuML.pyc config.pkl ensemble_model.pkl error tmp_scatter* normalize_dict.pkl
cp ./* ./Baseline
python runModel.py -l 1 -f general_gradients.cfg
./make_responses
./postproc
