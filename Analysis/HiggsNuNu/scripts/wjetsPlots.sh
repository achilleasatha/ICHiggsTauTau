#!/bin/sh
PARAMS=./scripts/Params.dat


###### n_muons
./bin/ControlPlots --cfg=scripts/controlPlot.cfg  \
    --plot_name="n_muons"  --x_axis_label="Number of muons" \
    --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=5 \
    --y_axis_min=0.01 --extra_pad=10000 \
    --rebin=1 \
    --norm_bins=false --verbose=false \
    --log_y=true \
    --blind=true --x_blind_min=1 --x_blind_max=0 \
    --paramfile=$PARAMS
./bin/ControlPlots --cfg=scripts/controlPlot.cfg  \
    --plot_name="n_muons"  --x_axis_label="Number of muons" \
    --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=5 \
    --y_axis_min=0.01 --extra_pad=2 \
    --rebin=1 \
    --norm_bins=false --verbose=false \
    --log_y=false \
    --blind=true --x_blind_min=1 --x_blind_max=0 \
    --paramfile=$PARAMS

###### mT
./bin/ControlPlots --cfg=scripts/controlPlot.cfg  \
    --plot_name="mt_munu"  --x_axis_label="m_{T} (GeV)" \
    --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=440 \
    --y_axis_min=0.01 --extra_pad=10000 \
    --rebin=20 \
    --norm_bins=false --verbose=false \
    --log_y=true \
    --blind=true --x_blind_min=1 --x_blind_max=0 \
    --paramfile=$PARAMS

###### mupt_1
./bin/ControlPlots --cfg=scripts/controlPlot.cfg  \
    --plot_name="mupt_1"  --x_axis_label="Leading muon p_{T} [GeV]" \
    --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=400 \
    --y_axis_min=0.01 --extra_pad=10000 \
    --rebin=40 \
    --norm_bins=false --verbose=false \
    --log_y=true \
    --blind=true --x_blind_min=1 --x_blind_max=0 \
    --paramfile=$PARAMS

###### mueta_1
./bin/ControlPlots --cfg=scripts/controlPlot.cfg  \
    --plot_name="mueta_1"  --x_axis_label="Leading muon #eta" \
    --custom_x_axis_range=true --x_axis_min=-3 --x_axis_max=3 \
    --y_axis_min=0.01 --extra_pad=10000 \
    --rebin=5 \
    --norm_bins=false --verbose=false \
    --log_y=true \
    --blind=true --x_blind_min=1 --x_blind_max=0 \
    --paramfile=$PARAMS

###### met
./bin/ControlPlots --cfg=scripts/controlPlot.cfg  \
    --plot_name="met_nomuons"  --x_axis_label="PF MET (GeV)" \
    --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=1000 \
    --y_axis_min=0.01 --extra_pad=10000 \
    --rebin=20 \
    --norm_bins=false \
    --log_y=true \
    --blind=true --x_blind_min=1 --x_blind_max=0 \
    --paramfile=$PARAMS
