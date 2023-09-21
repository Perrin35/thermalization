OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7744301) q[0];
sx q[0];
rz(-0.91355938) q[0];
sx q[0];
rz(1.4120742) q[0];
rz(-2.9867759) q[1];
sx q[1];
rz(-0.59564367) q[1];
sx q[1];
rz(-1.6593978) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9528113) q[0];
sx q[0];
rz(-2.6872098) q[0];
sx q[0];
rz(0.36034583) q[0];
rz(2.8843845) q[2];
sx q[2];
rz(-1.4423941) q[2];
sx q[2];
rz(0.53127015) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9285779) q[1];
sx q[1];
rz(-1.1482114) q[1];
sx q[1];
rz(-1.0282474) q[1];
rz(-pi) q[2];
rz(2.1336018) q[3];
sx q[3];
rz(-1.4853012) q[3];
sx q[3];
rz(0.29719719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.98510629) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(-2.2757754) q[2];
rz(0.95430294) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(1.2877864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99825478) q[0];
sx q[0];
rz(-1.7049494) q[0];
sx q[0];
rz(-0.026219333) q[0];
rz(-1.6014618) q[1];
sx q[1];
rz(-1.5988348) q[1];
sx q[1];
rz(2.1781133) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022097691) q[0];
sx q[0];
rz(-0.95675981) q[0];
sx q[0];
rz(-0.0028145785) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33032592) q[2];
sx q[2];
rz(-1.0268372) q[2];
sx q[2];
rz(-0.75737539) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3669489) q[1];
sx q[1];
rz(-1.0732713) q[1];
sx q[1];
rz(0.36555396) q[1];
x q[2];
rz(1.2198592) q[3];
sx q[3];
rz(-1.3388472) q[3];
sx q[3];
rz(1.9333145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6271237) q[2];
sx q[2];
rz(-1.1273948) q[2];
sx q[2];
rz(0.13452402) q[2];
rz(0.7450122) q[3];
sx q[3];
rz(-0.22694215) q[3];
sx q[3];
rz(0.9427332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9298252) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(2.3441558) q[0];
rz(1.047661) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(2.581596) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3968351) q[0];
sx q[0];
rz(-0.43733894) q[0];
sx q[0];
rz(2.0646981) q[0];
rz(-2.838344) q[2];
sx q[2];
rz(-1.5504642) q[2];
sx q[2];
rz(0.081239935) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9251688) q[1];
sx q[1];
rz(-2.7270626) q[1];
sx q[1];
rz(-0.51149909) q[1];
rz(-pi) q[2];
rz(-2.0542164) q[3];
sx q[3];
rz(-1.1421575) q[3];
sx q[3];
rz(-2.7408858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.75227633) q[2];
sx q[2];
rz(-1.9439149) q[2];
sx q[2];
rz(-0.17253549) q[2];
rz(0.98207384) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(-2.0836232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1383706) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(-0.28451434) q[0];
rz(-2.8248887) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(1.2987312) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0239149) q[0];
sx q[0];
rz(-1.0751343) q[0];
sx q[0];
rz(-0.25650521) q[0];
rz(-3.1403055) q[2];
sx q[2];
rz(-2.3372071) q[2];
sx q[2];
rz(0.019891642) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.270077) q[1];
sx q[1];
rz(-1.3844826) q[1];
sx q[1];
rz(2.1427076) q[1];
rz(-pi) q[2];
rz(-0.11233791) q[3];
sx q[3];
rz(-1.2445407) q[3];
sx q[3];
rz(2.5374967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5359042) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(0.38468012) q[2];
rz(2.3875333) q[3];
sx q[3];
rz(-2.0850756) q[3];
sx q[3];
rz(-1.4543021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5383179) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(1.3866562) q[0];
rz(-0.23100135) q[1];
sx q[1];
rz(-1.341154) q[1];
sx q[1];
rz(2.8447661) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41243991) q[0];
sx q[0];
rz(-1.531633) q[0];
sx q[0];
rz(2.5705283) q[0];
x q[1];
rz(2.9752521) q[2];
sx q[2];
rz(-0.74976774) q[2];
sx q[2];
rz(1.9094085) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0823114) q[1];
sx q[1];
rz(-2.0356405) q[1];
sx q[1];
rz(2.6962198) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0389901) q[3];
sx q[3];
rz(-0.71628621) q[3];
sx q[3];
rz(2.7308381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37830535) q[2];
sx q[2];
rz(-1.3102691) q[2];
sx q[2];
rz(-0.39247593) q[2];
rz(1.9893507) q[3];
sx q[3];
rz(-2.4270054) q[3];
sx q[3];
rz(0.31744441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0157938) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(-0.75138599) q[0];
rz(-1.3279351) q[1];
sx q[1];
rz(-1.8782047) q[1];
sx q[1];
rz(0.60633916) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7327001) q[0];
sx q[0];
rz(-0.048763976) q[0];
sx q[0];
rz(-0.34838895) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0688071) q[2];
sx q[2];
rz(-0.74337465) q[2];
sx q[2];
rz(2.5943287) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9669173) q[1];
sx q[1];
rz(-1.7511909) q[1];
sx q[1];
rz(2.128152) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.353225) q[3];
sx q[3];
rz(-1.4975582) q[3];
sx q[3];
rz(-2.6069802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5027344) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(-1.139337) q[2];
rz(-1.6566488) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(0.10425723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5320324) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(2.6334921) q[0];
rz(-1.5628901) q[1];
sx q[1];
rz(-1.0888313) q[1];
sx q[1];
rz(-2.3513444) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3944053) q[0];
sx q[0];
rz(-2.1571775) q[0];
sx q[0];
rz(1.9275083) q[0];
rz(0.21364084) q[2];
sx q[2];
rz(-1.5822615) q[2];
sx q[2];
rz(1.2917047) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7888469) q[1];
sx q[1];
rz(-0.48301304) q[1];
sx q[1];
rz(-0.61139815) q[1];
rz(1.1931476) q[3];
sx q[3];
rz(-2.1834063) q[3];
sx q[3];
rz(-2.2539504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.053085176) q[2];
sx q[2];
rz(-2.6997456) q[2];
sx q[2];
rz(-1.4132168) q[2];
rz(1.6648071) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(-0.061554519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4247894) q[0];
sx q[0];
rz(-3.1112818) q[0];
sx q[0];
rz(-2.0943663) q[0];
rz(2.5324902) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(-1.3887127) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0213288) q[0];
sx q[0];
rz(-1.5409894) q[0];
sx q[0];
rz(-2.1304312) q[0];
x q[1];
rz(-1.9954761) q[2];
sx q[2];
rz(-2.0600024) q[2];
sx q[2];
rz(0.81791544) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3988406) q[1];
sx q[1];
rz(-1.5821777) q[1];
sx q[1];
rz(3.0497754) q[1];
rz(-1.1561398) q[3];
sx q[3];
rz(-2.4932043) q[3];
sx q[3];
rz(0.63601953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9528815) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(0.88225538) q[2];
rz(-1.7404209) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(-1.2004948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27154487) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(-1.4260938) q[0];
rz(3.0601314) q[1];
sx q[1];
rz(-1.9790244) q[1];
sx q[1];
rz(0.55823278) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4687913) q[0];
sx q[0];
rz(-1.1835915) q[0];
sx q[0];
rz(-2.0331435) q[0];
x q[1];
rz(-0.054106601) q[2];
sx q[2];
rz(-1.6781085) q[2];
sx q[2];
rz(3.0215614) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8856814) q[1];
sx q[1];
rz(-2.3951888) q[1];
sx q[1];
rz(1.4766775) q[1];
x q[2];
rz(2.8006323) q[3];
sx q[3];
rz(-2.6302359) q[3];
sx q[3];
rz(-2.2380059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8490863) q[2];
sx q[2];
rz(-1.2693274) q[2];
sx q[2];
rz(-0.212184) q[2];
rz(0.21197453) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(1.9395444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6367209) q[0];
sx q[0];
rz(-0.8240521) q[0];
sx q[0];
rz(-1.5378392) q[0];
rz(0.82540712) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(2.6182981) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035318035) q[0];
sx q[0];
rz(-0.9629074) q[0];
sx q[0];
rz(1.4950698) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0977544) q[2];
sx q[2];
rz(-0.99928108) q[2];
sx q[2];
rz(0.29585719) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8955252) q[1];
sx q[1];
rz(-1.2125891) q[1];
sx q[1];
rz(-1.2342865) q[1];
rz(-pi) q[2];
rz(0.75533112) q[3];
sx q[3];
rz(-2.8108201) q[3];
sx q[3];
rz(-1.8473234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.837073) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(-2.8722897) q[2];
rz(-2.6473911) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(1.0860898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8158648) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(-1.4670463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(-1.5291443) q[2];
sx q[2];
rz(-0.26298444) q[2];
sx q[2];
rz(2.0900805) q[2];
rz(-0.79694637) q[3];
sx q[3];
rz(-1.7399825) q[3];
sx q[3];
rz(0.83132838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
