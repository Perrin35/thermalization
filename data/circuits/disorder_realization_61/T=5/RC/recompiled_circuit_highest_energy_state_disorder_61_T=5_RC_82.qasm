OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.8340303) q[0];
sx q[0];
rz(-2.7355255) q[0];
sx q[0];
rz(1.9598444) q[0];
rz(1.740068) q[1];
sx q[1];
rz(-0.47456804) q[1];
sx q[1];
rz(0.13303703) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1057464) q[0];
sx q[0];
rz(-1.3091358) q[0];
sx q[0];
rz(-2.2308971) q[0];
rz(0.51271397) q[2];
sx q[2];
rz(-1.9641335) q[2];
sx q[2];
rz(-2.8296997) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.052357167) q[1];
sx q[1];
rz(-1.3468139) q[1];
sx q[1];
rz(1.2473945) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0361133) q[3];
sx q[3];
rz(-2.4991192) q[3];
sx q[3];
rz(1.5758663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1349606) q[2];
sx q[2];
rz(-1.0198318) q[2];
sx q[2];
rz(2.4572065) q[2];
rz(1.1413752) q[3];
sx q[3];
rz(-2.8751774) q[3];
sx q[3];
rz(3.0548837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77966493) q[0];
sx q[0];
rz(-2.4517038) q[0];
sx q[0];
rz(2.6806114) q[0];
rz(-2.0224679) q[1];
sx q[1];
rz(-2.7179167) q[1];
sx q[1];
rz(2.8295595) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082816007) q[0];
sx q[0];
rz(-0.030555336) q[0];
sx q[0];
rz(-2.729134) q[0];
x q[1];
rz(-0.076032555) q[2];
sx q[2];
rz(-2.3302148) q[2];
sx q[2];
rz(1.3181612) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4496135) q[1];
sx q[1];
rz(-0.17573729) q[1];
sx q[1];
rz(-1.6040317) q[1];
rz(-pi) q[2];
rz(2.2226238) q[3];
sx q[3];
rz(-1.3462092) q[3];
sx q[3];
rz(-0.3995527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1877039) q[2];
sx q[2];
rz(-1.1053332) q[2];
sx q[2];
rz(0.40404955) q[2];
rz(0.36014253) q[3];
sx q[3];
rz(-1.6481684) q[3];
sx q[3];
rz(2.4586316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3563511) q[0];
sx q[0];
rz(-2.1389565) q[0];
sx q[0];
rz(0.32401618) q[0];
rz(-2.0815381) q[1];
sx q[1];
rz(-2.8476604) q[1];
sx q[1];
rz(0.70045984) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4840317) q[0];
sx q[0];
rz(-0.23398384) q[0];
sx q[0];
rz(-2.3768539) q[0];
rz(-pi) q[1];
rz(-1.9159928) q[2];
sx q[2];
rz(-1.3681865) q[2];
sx q[2];
rz(2.4576296) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6571275) q[1];
sx q[1];
rz(-1.336369) q[1];
sx q[1];
rz(-0.77204324) q[1];
rz(-pi) q[2];
rz(1.9869119) q[3];
sx q[3];
rz(-2.6857102) q[3];
sx q[3];
rz(2.1385764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6274274) q[2];
sx q[2];
rz(-0.24719396) q[2];
sx q[2];
rz(2.3449281) q[2];
rz(-1.1967777) q[3];
sx q[3];
rz(-1.5408885) q[3];
sx q[3];
rz(0.59320199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45348039) q[0];
sx q[0];
rz(-2.1402833) q[0];
sx q[0];
rz(-1.8185115) q[0];
rz(2.6755013) q[1];
sx q[1];
rz(-0.90704647) q[1];
sx q[1];
rz(1.1493433) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3484584) q[0];
sx q[0];
rz(-2.0597947) q[0];
sx q[0];
rz(-2.6455946) q[0];
rz(2.7485899) q[2];
sx q[2];
rz(-2.1519024) q[2];
sx q[2];
rz(2.4266134) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.016952) q[1];
sx q[1];
rz(-1.6158544) q[1];
sx q[1];
rz(2.9393815) q[1];
x q[2];
rz(1.1584267) q[3];
sx q[3];
rz(-2.7791808) q[3];
sx q[3];
rz(2.6965237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.65679067) q[2];
sx q[2];
rz(-0.66978407) q[2];
sx q[2];
rz(2.316324) q[2];
rz(-1.2800062) q[3];
sx q[3];
rz(-1.620159) q[3];
sx q[3];
rz(2.4204204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44329062) q[0];
sx q[0];
rz(-1.2606786) q[0];
sx q[0];
rz(-1.3662421) q[0];
rz(-2.0514936) q[1];
sx q[1];
rz(-1.6366199) q[1];
sx q[1];
rz(1.3390138) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96025673) q[0];
sx q[0];
rz(-0.68013817) q[0];
sx q[0];
rz(-0.28016092) q[0];
rz(-pi) q[1];
rz(2.2721185) q[2];
sx q[2];
rz(-0.83842248) q[2];
sx q[2];
rz(-2.7039607) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4046062) q[1];
sx q[1];
rz(-1.4530208) q[1];
sx q[1];
rz(-0.21356689) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3424266) q[3];
sx q[3];
rz(-0.46119565) q[3];
sx q[3];
rz(2.9950664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.94685405) q[2];
sx q[2];
rz(-1.7297435) q[2];
sx q[2];
rz(-0.33352387) q[2];
rz(-0.58878318) q[3];
sx q[3];
rz(-0.47313658) q[3];
sx q[3];
rz(0.80140448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3280585) q[0];
sx q[0];
rz(-1.3575587) q[0];
sx q[0];
rz(-1.300977) q[0];
rz(2.1133555) q[1];
sx q[1];
rz(-0.49003777) q[1];
sx q[1];
rz(-2.6992544) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16602941) q[0];
sx q[0];
rz(-1.1021758) q[0];
sx q[0];
rz(-1.2965078) q[0];
x q[1];
rz(-0.84965923) q[2];
sx q[2];
rz(-2.4772212) q[2];
sx q[2];
rz(-2.0400782) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.29242) q[1];
sx q[1];
rz(-1.994881) q[1];
sx q[1];
rz(1.7522274) q[1];
x q[2];
rz(-1.4668999) q[3];
sx q[3];
rz(-1.2524309) q[3];
sx q[3];
rz(-1.1863134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3698795) q[2];
sx q[2];
rz(-1.5481202) q[2];
sx q[2];
rz(-2.6698574) q[2];
rz(-1.4419904) q[3];
sx q[3];
rz(-2.580018) q[3];
sx q[3];
rz(-2.877318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0278552) q[0];
sx q[0];
rz(-1.4677445) q[0];
sx q[0];
rz(2.9748528) q[0];
rz(-2.3475504) q[1];
sx q[1];
rz(-1.9414976) q[1];
sx q[1];
rz(0.099743191) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31064311) q[0];
sx q[0];
rz(-1.2590053) q[0];
sx q[0];
rz(-0.98592178) q[0];
rz(-pi) q[1];
rz(-2.1879467) q[2];
sx q[2];
rz(-1.8815478) q[2];
sx q[2];
rz(1.0342395) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.21799) q[1];
sx q[1];
rz(-1.0031219) q[1];
sx q[1];
rz(-1.7204549) q[1];
rz(-pi) q[2];
rz(-2.2325956) q[3];
sx q[3];
rz(-1.2332877) q[3];
sx q[3];
rz(2.8765529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.86228592) q[2];
sx q[2];
rz(-2.0760832) q[2];
sx q[2];
rz(0.12969895) q[2];
rz(1.8779523) q[3];
sx q[3];
rz(-1.940515) q[3];
sx q[3];
rz(-0.14550801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.170914) q[0];
sx q[0];
rz(-2.2190974) q[0];
sx q[0];
rz(-0.42651919) q[0];
rz(2.049394) q[1];
sx q[1];
rz(-2.009232) q[1];
sx q[1];
rz(2.138864) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9364081) q[0];
sx q[0];
rz(-2.4800088) q[0];
sx q[0];
rz(-1.1422045) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4187047) q[2];
sx q[2];
rz(-1.6405511) q[2];
sx q[2];
rz(-0.43089128) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.59834023) q[1];
sx q[1];
rz(-1.9045254) q[1];
sx q[1];
rz(-2.0327492) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9055481) q[3];
sx q[3];
rz(-2.0147284) q[3];
sx q[3];
rz(-0.027396552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2670474) q[2];
sx q[2];
rz(-2.5858877) q[2];
sx q[2];
rz(0.12346539) q[2];
rz(0.63001436) q[3];
sx q[3];
rz(-1.29653) q[3];
sx q[3];
rz(-0.71648487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89710871) q[0];
sx q[0];
rz(-2.4429595) q[0];
sx q[0];
rz(-0.83936349) q[0];
rz(-0.15946236) q[1];
sx q[1];
rz(-0.99226743) q[1];
sx q[1];
rz(0.92963591) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2004505) q[0];
sx q[0];
rz(-1.8243891) q[0];
sx q[0];
rz(0.7483878) q[0];
rz(1.4484181) q[2];
sx q[2];
rz(-2.158383) q[2];
sx q[2];
rz(-1.1002968) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8690455) q[1];
sx q[1];
rz(-2.2490977) q[1];
sx q[1];
rz(-0.1609623) q[1];
x q[2];
rz(0.93299039) q[3];
sx q[3];
rz(-1.7562885) q[3];
sx q[3];
rz(2.958667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.9440426) q[2];
sx q[2];
rz(-1.9746747) q[2];
sx q[2];
rz(-2.3738532) q[2];
rz(-0.36457101) q[3];
sx q[3];
rz(-1.3837827) q[3];
sx q[3];
rz(3.072123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0070873) q[0];
sx q[0];
rz(-2.5275079) q[0];
sx q[0];
rz(-0.8836723) q[0];
rz(-0.20835749) q[1];
sx q[1];
rz(-2.6388984) q[1];
sx q[1];
rz(1.3349894) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17598957) q[0];
sx q[0];
rz(-1.5056247) q[0];
sx q[0];
rz(1.6287281) q[0];
rz(2.8284188) q[2];
sx q[2];
rz(-0.54944456) q[2];
sx q[2];
rz(-0.012481364) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6328278) q[1];
sx q[1];
rz(-2.0109051) q[1];
sx q[1];
rz(0.97036636) q[1];
rz(-2.317071) q[3];
sx q[3];
rz(-1.0484277) q[3];
sx q[3];
rz(2.8461412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16365446) q[2];
sx q[2];
rz(-1.6341219) q[2];
sx q[2];
rz(0.046023544) q[2];
rz(2.9764791) q[3];
sx q[3];
rz(-0.29296994) q[3];
sx q[3];
rz(1.2142115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56275) q[0];
sx q[0];
rz(-3.0414707) q[0];
sx q[0];
rz(1.1895251) q[0];
rz(0.1651172) q[1];
sx q[1];
rz(-1.6830291) q[1];
sx q[1];
rz(2.8370636) q[1];
rz(-0.61928673) q[2];
sx q[2];
rz(-2.613667) q[2];
sx q[2];
rz(2.1783301) q[2];
rz(2.3143309) q[3];
sx q[3];
rz(-0.65613366) q[3];
sx q[3];
rz(-0.14449672) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
