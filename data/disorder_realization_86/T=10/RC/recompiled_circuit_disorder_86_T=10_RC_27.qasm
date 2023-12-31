OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3671626) q[0];
sx q[0];
rz(-2.2280333) q[0];
sx q[0];
rz(1.7295184) q[0];
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
rz(2.1971561) q[0];
sx q[0];
rz(-1.4154139) q[0];
sx q[0];
rz(0.42874254) q[0];
x q[1];
rz(2.8843845) q[2];
sx q[2];
rz(-1.4423941) q[2];
sx q[2];
rz(0.53127015) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.60018051) q[1];
sx q[1];
rz(-2.0611079) q[1];
sx q[1];
rz(2.6580826) q[1];
rz(-pi) q[2];
rz(3.0406038) q[3];
sx q[3];
rz(-2.1312993) q[3];
sx q[3];
rz(-1.8141754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1564864) q[2];
sx q[2];
rz(-2.6323695) q[2];
sx q[2];
rz(0.86581725) q[2];
rz(-0.95430294) q[3];
sx q[3];
rz(-1.538397) q[3];
sx q[3];
rz(1.2877864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1433379) q[0];
sx q[0];
rz(-1.7049494) q[0];
sx q[0];
rz(-0.026219333) q[0];
rz(1.5401309) q[1];
sx q[1];
rz(-1.5988348) q[1];
sx q[1];
rz(-0.96347934) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5945157) q[0];
sx q[0];
rz(-1.5684959) q[0];
sx q[0];
rz(2.1848347) q[0];
rz(-pi) q[1];
rz(-2.1396779) q[2];
sx q[2];
rz(-1.2895673) q[2];
sx q[2];
rz(-0.98904726) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7746437) q[1];
sx q[1];
rz(-1.0732713) q[1];
sx q[1];
rz(0.36555396) q[1];
x q[2];
rz(-1.2198592) q[3];
sx q[3];
rz(-1.8027455) q[3];
sx q[3];
rz(1.9333145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5144689) q[2];
sx q[2];
rz(-1.1273948) q[2];
sx q[2];
rz(3.0070686) q[2];
rz(0.7450122) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(-0.9427332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298252) q[0];
sx q[0];
rz(-0.38914248) q[0];
sx q[0];
rz(2.3441558) q[0];
rz(-2.0939317) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(-0.55999666) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3968351) q[0];
sx q[0];
rz(-2.7042537) q[0];
sx q[0];
rz(-2.0646981) q[0];
rz(-pi) q[1];
rz(0.067990818) q[2];
sx q[2];
rz(-2.8376841) q[2];
sx q[2];
rz(-1.7169203) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8289889) q[1];
sx q[1];
rz(-1.3723515) q[1];
sx q[1];
rz(-0.36638422) q[1];
rz(2.3476944) q[3];
sx q[3];
rz(-2.5069935) q[3];
sx q[3];
rz(-1.8397699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.75227633) q[2];
sx q[2];
rz(-1.1976778) q[2];
sx q[2];
rz(-0.17253549) q[2];
rz(0.98207384) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(1.0579695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1383706) q[0];
sx q[0];
rz(-2.0439742) q[0];
sx q[0];
rz(-0.28451434) q[0];
rz(0.31670397) q[1];
sx q[1];
rz(-2.7088294) q[1];
sx q[1];
rz(1.8428615) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0239149) q[0];
sx q[0];
rz(-1.0751343) q[0];
sx q[0];
rz(-2.8850874) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1403055) q[2];
sx q[2];
rz(-2.3372071) q[2];
sx q[2];
rz(3.121701) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5609834) q[1];
sx q[1];
rz(-2.1316075) q[1];
sx q[1];
rz(2.9210655) q[1];
rz(-pi) q[2];
rz(1.2426162) q[3];
sx q[3];
rz(-1.4644074) q[3];
sx q[3];
rz(-1.002841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6056885) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(-2.7569125) q[2];
rz(-2.3875333) q[3];
sx q[3];
rz(-2.0850756) q[3];
sx q[3];
rz(-1.6872905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5383179) q[0];
sx q[0];
rz(-2.1544927) q[0];
sx q[0];
rz(1.3866562) q[0];
rz(-2.9105913) q[1];
sx q[1];
rz(-1.341154) q[1];
sx q[1];
rz(-2.8447661) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2192229) q[0];
sx q[0];
rz(-0.57225675) q[0];
sx q[0];
rz(0.072364307) q[0];
x q[1];
rz(-2.3987531) q[2];
sx q[2];
rz(-1.683871) q[2];
sx q[2];
rz(-2.6807221) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0592812) q[1];
sx q[1];
rz(-2.0356405) q[1];
sx q[1];
rz(-2.6962198) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71367587) q[3];
sx q[3];
rz(-1.6380966) q[3];
sx q[3];
rz(-1.0825368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.37830535) q[2];
sx q[2];
rz(-1.3102691) q[2];
sx q[2];
rz(-2.7491167) q[2];
rz(-1.9893507) q[3];
sx q[3];
rz(-2.4270054) q[3];
sx q[3];
rz(2.8241482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1257989) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(-0.75138599) q[0];
rz(1.8136576) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(2.5352535) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0814708) q[0];
sx q[0];
rz(-1.524964) q[0];
sx q[0];
rz(1.5541374) q[0];
rz(-pi) q[1];
rz(2.7251284) q[2];
sx q[2];
rz(-0.93566862) q[2];
sx q[2];
rz(-0.093402775) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11583466) q[1];
sx q[1];
rz(-2.5587213) q[1];
sx q[1];
rz(1.238766) q[1];
x q[2];
rz(-3.0665928) q[3];
sx q[3];
rz(-1.787775) q[3];
sx q[3];
rz(-1.0523588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.63885826) q[2];
sx q[2];
rz(-1.0417465) q[2];
sx q[2];
rz(2.0022557) q[2];
rz(1.4849439) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(-3.0373354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5320324) q[0];
sx q[0];
rz(-0.74815265) q[0];
sx q[0];
rz(0.50810057) q[0];
rz(1.5628901) q[1];
sx q[1];
rz(-1.0888313) q[1];
sx q[1];
rz(-0.79024822) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3944053) q[0];
sx q[0];
rz(-0.98441511) q[0];
sx q[0];
rz(1.2140843) q[0];
x q[1];
rz(1.5825282) q[2];
sx q[2];
rz(-1.3571697) q[2];
sx q[2];
rz(2.8649883) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1192757) q[1];
sx q[1];
rz(-1.1806618) q[1];
sx q[1];
rz(-1.2783865) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4942057) q[3];
sx q[3];
rz(-1.8772519) q[3];
sx q[3];
rz(2.2341773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.053085176) q[2];
sx q[2];
rz(-2.6997456) q[2];
sx q[2];
rz(-1.7283758) q[2];
rz(1.6648071) q[3];
sx q[3];
rz(-2.1063185) q[3];
sx q[3];
rz(0.061554519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7168032) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(2.0943663) q[0];
rz(0.60910243) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(-1.75288) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6435218) q[0];
sx q[0];
rz(-0.56034351) q[0];
sx q[0];
rz(1.5146921) q[0];
rz(-pi) q[1];
rz(2.6128204) q[2];
sx q[2];
rz(-1.1985589) q[2];
sx q[2];
rz(-0.96226529) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3988406) q[1];
sx q[1];
rz(-1.5594149) q[1];
sx q[1];
rz(-3.0497754) q[1];
x q[2];
rz(-0.96447585) q[3];
sx q[3];
rz(-1.8165605) q[3];
sx q[3];
rz(-0.59734674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1887112) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(2.2593373) q[2];
rz(-1.7404209) q[3];
sx q[3];
rz(-1.975235) q[3];
sx q[3];
rz(-1.9410979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8700478) q[0];
sx q[0];
rz(-2.7247868) q[0];
sx q[0];
rz(-1.7154988) q[0];
rz(3.0601314) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(-0.55823278) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24960625) q[0];
sx q[0];
rz(-2.5476646) q[0];
sx q[0];
rz(-0.83017613) q[0];
rz(-pi) q[1];
rz(-1.4633281) q[2];
sx q[2];
rz(-1.5170013) q[2];
sx q[2];
rz(-1.4449643) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0135632) q[1];
sx q[1];
rz(-2.3131144) q[1];
sx q[1];
rz(-0.086704266) q[1];
rz(-pi) q[2];
rz(0.48645143) q[3];
sx q[3];
rz(-1.7351741) q[3];
sx q[3];
rz(0.96729507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2925064) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(2.9294087) q[2];
rz(2.9296181) q[3];
sx q[3];
rz(-0.68325716) q[3];
sx q[3];
rz(-1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50487173) q[0];
sx q[0];
rz(-0.8240521) q[0];
sx q[0];
rz(1.5378392) q[0];
rz(2.3161855) q[1];
sx q[1];
rz(-2.4688265) q[1];
sx q[1];
rz(2.6182981) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4921724) q[0];
sx q[0];
rz(-1.6329375) q[0];
sx q[0];
rz(-2.5323575) q[0];
rz(1.5027572) q[2];
sx q[2];
rz(-2.5685852) q[2];
sx q[2];
rz(-0.21493658) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.20269468) q[1];
sx q[1];
rz(-1.885186) q[1];
sx q[1];
rz(0.37757341) q[1];
rz(-pi) q[2];
rz(2.3862615) q[3];
sx q[3];
rz(-0.33077251) q[3];
sx q[3];
rz(1.2942693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3045197) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(2.8722897) q[2];
rz(-2.6473911) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(1.0860898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8158648) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(1.4670463) q[1];
sx q[1];
rz(-2.84927) q[1];
sx q[1];
rz(-1.897859) q[1];
rz(1.6124484) q[2];
sx q[2];
rz(-0.26298444) q[2];
sx q[2];
rz(2.0900805) q[2];
rz(0.23444093) q[3];
sx q[3];
rz(-0.81080484) q[3];
sx q[3];
rz(-0.57639359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
