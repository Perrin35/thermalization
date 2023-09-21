OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6498123) q[0];
sx q[0];
rz(-0.28591135) q[0];
sx q[0];
rz(-2.6262992) q[0];
rz(4.4858785) q[1];
sx q[1];
rz(2.9872515) q[1];
sx q[1];
rz(6.8607688) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87969765) q[0];
sx q[0];
rz(-1.7356153) q[0];
sx q[0];
rz(-0.66189712) q[0];
rz(-pi) q[1];
rz(2.3876786) q[2];
sx q[2];
rz(-2.4587817) q[2];
sx q[2];
rz(-2.4726601) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2540993) q[1];
sx q[1];
rz(-1.3026397) q[1];
sx q[1];
rz(1.9894132) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64896119) q[3];
sx q[3];
rz(-1.0421841) q[3];
sx q[3];
rz(1.987207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.704533) q[2];
sx q[2];
rz(-1.5455064) q[2];
sx q[2];
rz(0.68721592) q[2];
rz(2.1263188) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(3.0190873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17094831) q[0];
sx q[0];
rz(-1.0785372) q[0];
sx q[0];
rz(1.8815536) q[0];
rz(-2.1353703) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(-2.2959183) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0802404) q[0];
sx q[0];
rz(-2.4363359) q[0];
sx q[0];
rz(2.4387226) q[0];
x q[1];
rz(0.13375608) q[2];
sx q[2];
rz(-1.6998561) q[2];
sx q[2];
rz(1.2226906) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0452803) q[1];
sx q[1];
rz(-1.7573866) q[1];
sx q[1];
rz(0.50218302) q[1];
rz(-pi) q[2];
rz(0.50176974) q[3];
sx q[3];
rz(-1.0106777) q[3];
sx q[3];
rz(0.027651699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.78850293) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(2.6611924) q[2];
rz(-1.7885615) q[3];
sx q[3];
rz(-1.055911) q[3];
sx q[3];
rz(1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1903494) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(2.3685266) q[0];
rz(3.0103325) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(-2.0551596) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0057384) q[0];
sx q[0];
rz(-1.5691301) q[0];
sx q[0];
rz(2.0322331) q[0];
rz(-pi) q[1];
rz(-1.418872) q[2];
sx q[2];
rz(-0.65290367) q[2];
sx q[2];
rz(2.8002847) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.74124386) q[1];
sx q[1];
rz(-1.8596008) q[1];
sx q[1];
rz(-1.5653533) q[1];
x q[2];
rz(1.9589892) q[3];
sx q[3];
rz(-2.7079294) q[3];
sx q[3];
rz(0.38503034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1290258) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(-1.770299) q[2];
rz(2.7584372) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(-0.80254054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4521769) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(3.0932328) q[0];
rz(-2.9776749) q[1];
sx q[1];
rz(-0.36968958) q[1];
sx q[1];
rz(1.4455459) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4500344) q[0];
sx q[0];
rz(-1.3908747) q[0];
sx q[0];
rz(-0.66288373) q[0];
rz(-pi) q[1];
rz(0.38979001) q[2];
sx q[2];
rz(-2.9082657) q[2];
sx q[2];
rz(2.0248272) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85019894) q[1];
sx q[1];
rz(-1.351601) q[1];
sx q[1];
rz(1.6367957) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0936071) q[3];
sx q[3];
rz(-2.5683937) q[3];
sx q[3];
rz(1.5447865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0754898) q[2];
sx q[2];
rz(-1.3572071) q[2];
sx q[2];
rz(-1.0243105) q[2];
rz(1.5284437) q[3];
sx q[3];
rz(-1.5214835) q[3];
sx q[3];
rz(2.9083692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4399453) q[0];
sx q[0];
rz(-2.3174536) q[0];
sx q[0];
rz(-1.8540927) q[0];
rz(0.31907407) q[1];
sx q[1];
rz(-1.5417475) q[1];
sx q[1];
rz(-2.2873926) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6720649) q[0];
sx q[0];
rz(-0.85692353) q[0];
sx q[0];
rz(-3.1307427) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87128432) q[2];
sx q[2];
rz(-1.9565906) q[2];
sx q[2];
rz(-2.3988349) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.90480587) q[1];
sx q[1];
rz(-0.85610897) q[1];
sx q[1];
rz(0.76101117) q[1];
rz(-0.55027996) q[3];
sx q[3];
rz(-2.0645421) q[3];
sx q[3];
rz(2.9328336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95191082) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(2.5642776) q[2];
rz(-2.632085) q[3];
sx q[3];
rz(-2.7039492) q[3];
sx q[3];
rz(-0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0034870738) q[0];
sx q[0];
rz(-2.0697937) q[0];
sx q[0];
rz(3.0694718) q[0];
rz(2.0347319) q[1];
sx q[1];
rz(-0.51263428) q[1];
sx q[1];
rz(-3.0153826) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.135658) q[0];
sx q[0];
rz(-1.6023484) q[0];
sx q[0];
rz(-2.9289782) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1926786) q[2];
sx q[2];
rz(-0.83231229) q[2];
sx q[2];
rz(-0.2371012) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8773552) q[1];
sx q[1];
rz(-2.4009631) q[1];
sx q[1];
rz(-1.2901558) q[1];
rz(1.24228) q[3];
sx q[3];
rz(-1.7528755) q[3];
sx q[3];
rz(-1.3086705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.90298992) q[2];
sx q[2];
rz(-2.949252) q[2];
sx q[2];
rz(-2.3664756) q[2];
rz(2.3136247) q[3];
sx q[3];
rz(-0.29100806) q[3];
sx q[3];
rz(-1.1221788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5489952) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(3.0431842) q[0];
rz(-1.1920284) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(-0.55955204) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0347621) q[0];
sx q[0];
rz(-2.1666399) q[0];
sx q[0];
rz(-1.725127) q[0];
rz(-1.5188811) q[2];
sx q[2];
rz(-1.8318818) q[2];
sx q[2];
rz(-1.3494929) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2110062) q[1];
sx q[1];
rz(-1.483327) q[1];
sx q[1];
rz(-1.4409815) q[1];
rz(0.26569326) q[3];
sx q[3];
rz(-0.64722792) q[3];
sx q[3];
rz(-0.94097394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.45903912) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(2.7977978) q[2];
rz(-2.5750459) q[3];
sx q[3];
rz(-2.6930801) q[3];
sx q[3];
rz(-2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7664117) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(-2.4108316) q[0];
rz(-2.9991951) q[1];
sx q[1];
rz(-1.8715033) q[1];
sx q[1];
rz(0.87160814) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4260092) q[0];
sx q[0];
rz(-1.4953488) q[0];
sx q[0];
rz(1.5796608) q[0];
rz(-0.40558221) q[2];
sx q[2];
rz(-1.431258) q[2];
sx q[2];
rz(-1.9987193) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.826088) q[1];
sx q[1];
rz(-1.2102038) q[1];
sx q[1];
rz(0.70198595) q[1];
rz(-pi) q[2];
rz(-1.7730764) q[3];
sx q[3];
rz(-0.61198046) q[3];
sx q[3];
rz(-0.66093854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4153851) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(2.0020206) q[2];
rz(-1.4987882) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(-0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9386439) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(1.9198445) q[0];
rz(-2.9755759) q[1];
sx q[1];
rz(-1.8211726) q[1];
sx q[1];
rz(-1.5244012) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6451384) q[0];
sx q[0];
rz(-1.5418515) q[0];
sx q[0];
rz(1.4883947) q[0];
rz(-pi) q[1];
rz(1.7622856) q[2];
sx q[2];
rz(-2.0404437) q[2];
sx q[2];
rz(0.42275235) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7414654) q[1];
sx q[1];
rz(-1.2418081) q[1];
sx q[1];
rz(-0.16727438) q[1];
x q[2];
rz(-2.5025326) q[3];
sx q[3];
rz(-1.3165054) q[3];
sx q[3];
rz(-1.1921079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1200072) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(2.7900556) q[2];
rz(2.0848138) q[3];
sx q[3];
rz(-2.6119699) q[3];
sx q[3];
rz(-0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4979424) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(1.836401) q[0];
rz(-2.7611043) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(-0.25340733) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8587592) q[0];
sx q[0];
rz(-1.4677605) q[0];
sx q[0];
rz(-1.9673002) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2316522) q[2];
sx q[2];
rz(-1.7288496) q[2];
sx q[2];
rz(-0.16976419) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0903783) q[1];
sx q[1];
rz(-1.9131294) q[1];
sx q[1];
rz(1.315209) q[1];
rz(-pi) q[2];
rz(2.5640423) q[3];
sx q[3];
rz(-0.69637075) q[3];
sx q[3];
rz(-1.3814572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5499251) q[2];
sx q[2];
rz(-0.88576907) q[2];
sx q[2];
rz(-0.5029451) q[2];
rz(-2.2425966) q[3];
sx q[3];
rz(-1.2939724) q[3];
sx q[3];
rz(-1.9780654) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1702561) q[0];
sx q[0];
rz(-1.6032871) q[0];
sx q[0];
rz(0.26300318) q[0];
rz(-2.4304216) q[1];
sx q[1];
rz(-2.053459) q[1];
sx q[1];
rz(-1.4278535) q[1];
rz(-1.8297557) q[2];
sx q[2];
rz(-1.8428409) q[2];
sx q[2];
rz(0.98696282) q[2];
rz(-2.3861804) q[3];
sx q[3];
rz(-1.0676386) q[3];
sx q[3];
rz(-0.044015351) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];