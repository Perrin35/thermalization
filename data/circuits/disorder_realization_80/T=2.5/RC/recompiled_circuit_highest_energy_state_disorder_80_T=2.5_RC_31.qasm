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
rz(0.12240527) q[0];
sx q[0];
rz(-0.91949099) q[0];
sx q[0];
rz(1.9424196) q[0];
rz(0.1872669) q[1];
sx q[1];
rz(6.8254677) q[1];
sx q[1];
rz(10.944278) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8336415) q[0];
sx q[0];
rz(-2.184377) q[0];
sx q[0];
rz(1.475901) q[0];
x q[1];
rz(-1.1178944) q[2];
sx q[2];
rz(-1.9097966) q[2];
sx q[2];
rz(-1.1690559) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6295084) q[1];
sx q[1];
rz(-2.0633882) q[1];
sx q[1];
rz(-1.1118481) q[1];
rz(-1.5222129) q[3];
sx q[3];
rz(-1.3268927) q[3];
sx q[3];
rz(-0.78256202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2970994) q[2];
sx q[2];
rz(-0.62701925) q[2];
sx q[2];
rz(1.7681047) q[2];
rz(-2.3376236) q[3];
sx q[3];
rz(-1.6425902) q[3];
sx q[3];
rz(1.1871626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3961769) q[0];
sx q[0];
rz(-0.71393037) q[0];
sx q[0];
rz(0.3748689) q[0];
rz(0.51775852) q[1];
sx q[1];
rz(-1.9381783) q[1];
sx q[1];
rz(-1.7727324) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8086209) q[0];
sx q[0];
rz(-3.140641) q[0];
sx q[0];
rz(3.0406221) q[0];
x q[1];
rz(0.19646074) q[2];
sx q[2];
rz(-1.5443373) q[2];
sx q[2];
rz(-0.65633869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6102099) q[1];
sx q[1];
rz(-2.068559) q[1];
sx q[1];
rz(2.0396292) q[1];
rz(-pi) q[2];
rz(0.13624713) q[3];
sx q[3];
rz(-1.3582503) q[3];
sx q[3];
rz(0.023879026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0181197) q[2];
sx q[2];
rz(-0.70088434) q[2];
sx q[2];
rz(-2.4269721) q[2];
rz(-2.9361652) q[3];
sx q[3];
rz(-1.3222062) q[3];
sx q[3];
rz(1.8429168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4648723) q[0];
sx q[0];
rz(-2.1045852) q[0];
sx q[0];
rz(-0.49764693) q[0];
rz(2.5045555) q[1];
sx q[1];
rz(-2.5126845) q[1];
sx q[1];
rz(0.10674891) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7926559) q[0];
sx q[0];
rz(-1.35864) q[0];
sx q[0];
rz(-1.8769708) q[0];
rz(-pi) q[1];
rz(1.5791513) q[2];
sx q[2];
rz(-0.66973493) q[2];
sx q[2];
rz(-1.005203) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5046204) q[1];
sx q[1];
rz(-0.85831888) q[1];
sx q[1];
rz(-3.0709355) q[1];
rz(2.51744) q[3];
sx q[3];
rz(-2.660203) q[3];
sx q[3];
rz(0.91048756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.60483021) q[2];
sx q[2];
rz(-2.3894775) q[2];
sx q[2];
rz(-2.9467648) q[2];
rz(0.005006494) q[3];
sx q[3];
rz(-0.61401335) q[3];
sx q[3];
rz(2.3495638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1057338) q[0];
sx q[0];
rz(-2.6073313) q[0];
sx q[0];
rz(2.1015097) q[0];
rz(-2.9478574) q[1];
sx q[1];
rz(-1.6400784) q[1];
sx q[1];
rz(2.3949882) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9039041) q[0];
sx q[0];
rz(-2.6857566) q[0];
sx q[0];
rz(1.8346915) q[0];
rz(-pi) q[1];
rz(2.9992013) q[2];
sx q[2];
rz(-2.2658544) q[2];
sx q[2];
rz(2.5945239) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6317217) q[1];
sx q[1];
rz(-1.5190017) q[1];
sx q[1];
rz(-0.69433166) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6799404) q[3];
sx q[3];
rz(-1.8092844) q[3];
sx q[3];
rz(1.3772688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1497583) q[2];
sx q[2];
rz(-1.6712302) q[2];
sx q[2];
rz(0.20450083) q[2];
rz(1.9483942) q[3];
sx q[3];
rz(-0.85719332) q[3];
sx q[3];
rz(-0.96113718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1292773) q[0];
sx q[0];
rz(-0.17450541) q[0];
sx q[0];
rz(-2.0611064) q[0];
rz(-2.7642545) q[1];
sx q[1];
rz(-1.6308547) q[1];
sx q[1];
rz(-2.2468755) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5114488) q[0];
sx q[0];
rz(-2.7155657) q[0];
sx q[0];
rz(2.565284) q[0];
x q[1];
rz(2.3482125) q[2];
sx q[2];
rz(-2.0700244) q[2];
sx q[2];
rz(2.3069059) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1929568) q[1];
sx q[1];
rz(-2.0564449) q[1];
sx q[1];
rz(-2.5705277) q[1];
x q[2];
rz(-2.0438439) q[3];
sx q[3];
rz(-2.4425207) q[3];
sx q[3];
rz(-3.0300273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6696024) q[2];
sx q[2];
rz(-0.44595465) q[2];
sx q[2];
rz(3.0338244) q[2];
rz(-2.9660411) q[3];
sx q[3];
rz(-1.8864417) q[3];
sx q[3];
rz(-2.5811783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.046722978) q[0];
sx q[0];
rz(-0.51114285) q[0];
sx q[0];
rz(1.0119337) q[0];
rz(-1.5049505) q[1];
sx q[1];
rz(-2.4220146) q[1];
sx q[1];
rz(-2.124427) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1873774) q[0];
sx q[0];
rz(-0.34437505) q[0];
sx q[0];
rz(-2.1013174) q[0];
rz(-pi) q[1];
rz(-0.76061298) q[2];
sx q[2];
rz(-1.9991864) q[2];
sx q[2];
rz(1.9477159) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.575841) q[1];
sx q[1];
rz(-1.3804129) q[1];
sx q[1];
rz(-1.562005) q[1];
rz(-2.9842942) q[3];
sx q[3];
rz(-1.1264115) q[3];
sx q[3];
rz(-0.97287662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.73198685) q[2];
sx q[2];
rz(-0.66183949) q[2];
sx q[2];
rz(2.1448263) q[2];
rz(3.0766727) q[3];
sx q[3];
rz(-1.4109979) q[3];
sx q[3];
rz(-1.1499278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053059269) q[0];
sx q[0];
rz(-0.94045883) q[0];
sx q[0];
rz(0.9084107) q[0];
rz(2.9216595) q[1];
sx q[1];
rz(-1.5989774) q[1];
sx q[1];
rz(-1.251108) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.661099) q[0];
sx q[0];
rz(-2.0401883) q[0];
sx q[0];
rz(-0.29222699) q[0];
rz(-1.471453) q[2];
sx q[2];
rz(-1.2867498) q[2];
sx q[2];
rz(-0.63114511) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0323769) q[1];
sx q[1];
rz(-0.93689474) q[1];
sx q[1];
rz(-2.3486167) q[1];
x q[2];
rz(3.0165169) q[3];
sx q[3];
rz(-0.8440869) q[3];
sx q[3];
rz(-0.39966003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1088341) q[2];
sx q[2];
rz(-0.83403504) q[2];
sx q[2];
rz(-1.533482) q[2];
rz(1.1533302) q[3];
sx q[3];
rz(-2.0861552) q[3];
sx q[3];
rz(-1.4920894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.15025) q[0];
sx q[0];
rz(-2.7582176) q[0];
sx q[0];
rz(-2.2008994) q[0];
rz(-3.0062145) q[1];
sx q[1];
rz(-1.4043413) q[1];
sx q[1];
rz(-2.1167596) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.230119) q[0];
sx q[0];
rz(-2.1699454) q[0];
sx q[0];
rz(2.7658303) q[0];
x q[1];
rz(2.5338855) q[2];
sx q[2];
rz(-1.6204699) q[2];
sx q[2];
rz(-2.7756754) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.26014806) q[1];
sx q[1];
rz(-1.3760412) q[1];
sx q[1];
rz(1.4347726) q[1];
rz(-pi) q[2];
x q[2];
rz(3.023324) q[3];
sx q[3];
rz(-0.41384816) q[3];
sx q[3];
rz(-2.5914206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.7950874) q[2];
sx q[2];
rz(-2.2515209) q[2];
sx q[2];
rz(-2.4915462) q[2];
rz(2.5551689) q[3];
sx q[3];
rz(-1.4805877) q[3];
sx q[3];
rz(2.39095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8745678) q[0];
sx q[0];
rz(-0.50849193) q[0];
sx q[0];
rz(-1.1453999) q[0];
rz(1.3151431) q[1];
sx q[1];
rz(-1.2329085) q[1];
sx q[1];
rz(2.4536536) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2719587) q[0];
sx q[0];
rz(-0.42068538) q[0];
sx q[0];
rz(2.0482333) q[0];
x q[1];
rz(-2.1678989) q[2];
sx q[2];
rz(-1.8988238) q[2];
sx q[2];
rz(2.6231678) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.44837828) q[1];
sx q[1];
rz(-0.82469392) q[1];
sx q[1];
rz(-1.34971) q[1];
rz(-2.0618556) q[3];
sx q[3];
rz(-1.907067) q[3];
sx q[3];
rz(-0.46592135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.410586) q[2];
sx q[2];
rz(-2.035391) q[2];
sx q[2];
rz(2.2753687) q[2];
rz(-0.90803641) q[3];
sx q[3];
rz(-0.76483813) q[3];
sx q[3];
rz(-1.0959371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9662125) q[0];
sx q[0];
rz(-1.1840273) q[0];
sx q[0];
rz(2.7259735) q[0];
rz(0.79187727) q[1];
sx q[1];
rz(-1.917058) q[1];
sx q[1];
rz(-0.22722879) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89657518) q[0];
sx q[0];
rz(-1.5868717) q[0];
sx q[0];
rz(1.5590369) q[0];
rz(-pi) q[1];
rz(-0.80663075) q[2];
sx q[2];
rz(-1.5147594) q[2];
sx q[2];
rz(0.12412589) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.38293441) q[1];
sx q[1];
rz(-2.1442752) q[1];
sx q[1];
rz(-1.8783147) q[1];
rz(-2.0437743) q[3];
sx q[3];
rz(-2.0379582) q[3];
sx q[3];
rz(-0.43097365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8772584) q[2];
sx q[2];
rz(-0.72089583) q[2];
sx q[2];
rz(-0.43992511) q[2];
rz(-2.544493) q[3];
sx q[3];
rz(-2.4720981) q[3];
sx q[3];
rz(1.7841024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6482342) q[0];
sx q[0];
rz(-1.5536722) q[0];
sx q[0];
rz(-1.8552725) q[0];
rz(-1.413912) q[1];
sx q[1];
rz(-1.021011) q[1];
sx q[1];
rz(0.11722142) q[1];
rz(1.6261423) q[2];
sx q[2];
rz(-1.1641486) q[2];
sx q[2];
rz(2.138242) q[2];
rz(0.24735484) q[3];
sx q[3];
rz(-0.83360278) q[3];
sx q[3];
rz(-3.0616888) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
