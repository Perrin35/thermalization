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
rz(-1.7295184) q[0];
rz(-2.9867759) q[1];
sx q[1];
rz(5.6875416) q[1];
sx q[1];
rz(7.7653801) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55573758) q[0];
sx q[0];
rz(-1.9940388) q[0];
sx q[0];
rz(-1.7413571) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7035159) q[2];
sx q[2];
rz(-1.3157529) q[2];
sx q[2];
rz(2.1357352) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9285779) q[1];
sx q[1];
rz(-1.1482114) q[1];
sx q[1];
rz(-2.1133452) q[1];
rz(2.1336018) q[3];
sx q[3];
rz(-1.4853012) q[3];
sx q[3];
rz(-2.8443955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.98510629) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(-2.2757754) q[2];
rz(0.95430294) q[3];
sx q[3];
rz(-1.538397) q[3];
sx q[3];
rz(-1.2877864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99825478) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(3.1153733) q[0];
rz(1.5401309) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(-2.1781133) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.547077) q[0];
sx q[0];
rz(-1.5684959) q[0];
sx q[0];
rz(0.95675795) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0630402) q[2];
sx q[2];
rz(-0.62765593) q[2];
sx q[2];
rz(-0.17222675) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1184247) q[1];
sx q[1];
rz(-1.8903362) q[1];
sx q[1];
rz(1.0440786) q[1];
rz(-pi) q[2];
rz(0.24641896) q[3];
sx q[3];
rz(-1.911947) q[3];
sx q[3];
rz(2.6951172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5144689) q[2];
sx q[2];
rz(-2.0141979) q[2];
sx q[2];
rz(3.0070686) q[2];
rz(-0.7450122) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(-2.1988595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2117675) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(-2.3441558) q[0];
rz(1.047661) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(2.581596) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6277498) q[0];
sx q[0];
rz(-1.3686413) q[0];
sx q[0];
rz(1.1802799) q[0];
rz(3.0736018) q[2];
sx q[2];
rz(-0.30390856) q[2];
sx q[2];
rz(-1.7169203) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9251688) q[1];
sx q[1];
rz(-2.7270626) q[1];
sx q[1];
rz(0.51149909) q[1];
rz(-pi) q[2];
rz(2.6651354) q[3];
sx q[3];
rz(-2.0072848) q[3];
sx q[3];
rz(-0.95526327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75227633) q[2];
sx q[2];
rz(-1.1976778) q[2];
sx q[2];
rz(2.9690572) q[2];
rz(-2.1595188) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(1.0579695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.003222) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(-0.28451434) q[0];
rz(-2.8248887) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(1.2987312) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38628681) q[0];
sx q[0];
rz(-0.55314976) q[0];
sx q[0];
rz(-1.132071) q[0];
rz(-pi) q[1];
rz(3.1403055) q[2];
sx q[2];
rz(-0.80438559) q[2];
sx q[2];
rz(-3.121701) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.87151566) q[1];
sx q[1];
rz(-1.3844826) q[1];
sx q[1];
rz(2.1427076) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11233791) q[3];
sx q[3];
rz(-1.897052) q[3];
sx q[3];
rz(0.60409594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6056885) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(-2.7569125) q[2];
rz(0.7540594) q[3];
sx q[3];
rz(-2.0850756) q[3];
sx q[3];
rz(1.4543021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5383179) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(1.3866562) q[0];
rz(0.23100135) q[1];
sx q[1];
rz(-1.341154) q[1];
sx q[1];
rz(-2.8447661) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2192229) q[0];
sx q[0];
rz(-2.5693359) q[0];
sx q[0];
rz(-0.072364307) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9752521) q[2];
sx q[2];
rz(-2.3918249) q[2];
sx q[2];
rz(-1.2321842) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2424803) q[1];
sx q[1];
rz(-0.63226262) q[1];
sx q[1];
rz(-0.86109032) q[1];
rz(-pi) q[2];
rz(0.10260251) q[3];
sx q[3];
rz(-0.71628621) q[3];
sx q[3];
rz(-0.41075452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.37830535) q[2];
sx q[2];
rz(-1.8313235) q[2];
sx q[2];
rz(2.7491167) q[2];
rz(-1.1522419) q[3];
sx q[3];
rz(-2.4270054) q[3];
sx q[3];
rz(-2.8241482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1257989) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(-2.3902067) q[0];
rz(1.8136576) q[1];
sx q[1];
rz(-1.8782047) q[1];
sx q[1];
rz(-2.5352535) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5099112) q[0];
sx q[0];
rz(-1.5541549) q[0];
sx q[0];
rz(-3.095754) q[0];
rz(-pi) q[1];
rz(-2.0727856) q[2];
sx q[2];
rz(-0.74337465) q[2];
sx q[2];
rz(2.5943287) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9669173) q[1];
sx q[1];
rz(-1.7511909) q[1];
sx q[1];
rz(1.0134407) q[1];
rz(-pi) q[2];
rz(1.243152) q[3];
sx q[3];
rz(-2.9122105) q[3];
sx q[3];
rz(-0.71654191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5027344) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(-2.0022557) q[2];
rz(-1.6566488) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(-3.0373354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6095603) q[0];
sx q[0];
rz(-0.74815265) q[0];
sx q[0];
rz(-2.6334921) q[0];
rz(-1.5787026) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(0.79024822) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3944053) q[0];
sx q[0];
rz(-0.98441511) q[0];
sx q[0];
rz(1.9275083) q[0];
rz(-pi) q[1];
rz(3.087567) q[2];
sx q[2];
rz(-0.21394357) q[2];
sx q[2];
rz(-0.33188785) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7888469) q[1];
sx q[1];
rz(-2.6585796) q[1];
sx q[1];
rz(-2.5301945) q[1];
rz(-2.6584133) q[3];
sx q[3];
rz(-0.70671591) q[3];
sx q[3];
rz(-0.28373517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.053085176) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(-1.4132168) q[2];
rz(1.4767856) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(0.061554519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-1.7168032) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(-2.0943663) q[0];
rz(-0.60910243) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(-1.3887127) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0213288) q[0];
sx q[0];
rz(-1.6006032) q[0];
sx q[0];
rz(-1.0111615) q[0];
rz(-pi) q[1];
rz(0.52877229) q[2];
sx q[2];
rz(-1.1985589) q[2];
sx q[2];
rz(0.96226529) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3125004) q[1];
sx q[1];
rz(-1.6626076) q[1];
sx q[1];
rz(-1.5593668) q[1];
x q[2];
rz(-1.9854529) q[3];
sx q[3];
rz(-2.4932043) q[3];
sx q[3];
rz(-0.63601953) q[3];
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
rz(1.7404209) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(-1.9410979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27154487) q[0];
sx q[0];
rz(-2.7247868) q[0];
sx q[0];
rz(-1.7154988) q[0];
rz(0.081461279) q[1];
sx q[1];
rz(-1.9790244) q[1];
sx q[1];
rz(-0.55823278) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0840069) q[0];
sx q[0];
rz(-1.1450197) q[0];
sx q[0];
rz(0.42752479) q[0];
rz(-3.0874861) q[2];
sx q[2];
rz(-1.6781085) q[2];
sx q[2];
rz(0.12003128) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7575175) q[1];
sx q[1];
rz(-1.5069403) q[1];
sx q[1];
rz(0.82660316) q[1];
rz(-pi) q[2];
rz(-2.8006323) q[3];
sx q[3];
rz(-0.51135671) q[3];
sx q[3];
rz(-2.2380059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2925064) q[2];
sx q[2];
rz(-1.2693274) q[2];
sx q[2];
rz(-0.212184) q[2];
rz(0.21197453) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(-1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50487173) q[0];
sx q[0];
rz(-0.8240521) q[0];
sx q[0];
rz(-1.5378392) q[0];
rz(-0.82540712) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(0.5232946) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16738811) q[0];
sx q[0];
rz(-0.61199576) q[0];
sx q[0];
rz(-0.1083072) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5027572) q[2];
sx q[2];
rz(-2.5685852) q[2];
sx q[2];
rz(-0.21493658) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.20269468) q[1];
sx q[1];
rz(-1.2564066) q[1];
sx q[1];
rz(2.7640192) q[1];
rz(-1.80199) q[3];
sx q[3];
rz(-1.3320859) q[3];
sx q[3];
rz(-2.0774487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.837073) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(-2.8722897) q[2];
rz(-0.4942016) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(2.0555029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3257278) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(-1.6745463) q[1];
sx q[1];
rz(-2.84927) q[1];
sx q[1];
rz(-1.897859) q[1];
rz(-1.3080296) q[2];
sx q[2];
rz(-1.5599712) q[2];
sx q[2];
rz(0.47906265) q[2];
rz(-1.3310824) q[3];
sx q[3];
rz(-2.3532383) q[3];
sx q[3];
rz(-0.91010703) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];