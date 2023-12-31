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
rz(0.15481678) q[1];
sx q[1];
rz(-2.545949) q[1];
sx q[1];
rz(1.6593978) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94443653) q[0];
sx q[0];
rz(-1.4154139) q[0];
sx q[0];
rz(2.7128501) q[0];
rz(-pi) q[1];
rz(-1.4380768) q[2];
sx q[2];
rz(-1.8258397) q[2];
sx q[2];
rz(-2.1357352) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60018051) q[1];
sx q[1];
rz(-2.0611079) q[1];
sx q[1];
rz(-2.6580826) q[1];
rz(-pi) q[2];
rz(-3.0406038) q[3];
sx q[3];
rz(-2.1312993) q[3];
sx q[3];
rz(1.8141754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.98510629) q[2];
sx q[2];
rz(-2.6323695) q[2];
sx q[2];
rz(0.86581725) q[2];
rz(0.95430294) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(1.2877864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1433379) q[0];
sx q[0];
rz(-1.7049494) q[0];
sx q[0];
rz(-0.026219333) q[0];
rz(-1.5401309) q[1];
sx q[1];
rz(-1.5988348) q[1];
sx q[1];
rz(-2.1781133) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5945157) q[0];
sx q[0];
rz(-1.5730968) q[0];
sx q[0];
rz(-0.95675795) q[0];
x q[1];
rz(-2.8112667) q[2];
sx q[2];
rz(-1.0268372) q[2];
sx q[2];
rz(0.75737539) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7746437) q[1];
sx q[1];
rz(-2.0683214) q[1];
sx q[1];
rz(2.7760387) q[1];
rz(1.9217334) q[3];
sx q[3];
rz(-1.3388472) q[3];
sx q[3];
rz(1.2082781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6271237) q[2];
sx q[2];
rz(-1.1273948) q[2];
sx q[2];
rz(3.0070686) q[2];
rz(2.3965805) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(-2.1988595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298252) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(-2.3441558) q[0];
rz(1.047661) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(-0.55999666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0021734) q[0];
sx q[0];
rz(-1.1886485) q[0];
sx q[0];
rz(-2.9234773) q[0];
rz(-pi) q[1];
rz(0.30324869) q[2];
sx q[2];
rz(-1.5504642) q[2];
sx q[2];
rz(0.081239935) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.33369505) q[1];
sx q[1];
rz(-1.2119319) q[1];
sx q[1];
rz(1.3586504) q[1];
x q[2];
rz(-2.0542164) q[3];
sx q[3];
rz(-1.9994352) q[3];
sx q[3];
rz(-0.40070686) q[3];
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
rz(-0.98207384) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(-1.0579695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.003222) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(-2.8570783) q[0];
rz(-2.8248887) q[1];
sx q[1];
rz(-2.7088294) q[1];
sx q[1];
rz(1.8428615) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5772229) q[0];
sx q[0];
rz(-1.7958613) q[0];
sx q[0];
rz(1.0610915) q[0];
x q[1];
rz(-3.1403055) q[2];
sx q[2];
rz(-0.80438559) q[2];
sx q[2];
rz(-0.019891642) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1620996) q[1];
sx q[1];
rz(-0.59826189) q[1];
sx q[1];
rz(1.23566) q[1];
rz(3.0292547) q[3];
sx q[3];
rz(-1.2445407) q[3];
sx q[3];
rz(-0.60409594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6056885) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(2.7569125) q[2];
rz(-0.7540594) q[3];
sx q[3];
rz(-2.0850756) q[3];
sx q[3];
rz(1.6872905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5383179) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(1.7549365) q[0];
rz(2.9105913) q[1];
sx q[1];
rz(-1.341154) q[1];
sx q[1];
rz(2.8447661) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2192229) q[0];
sx q[0];
rz(-2.5693359) q[0];
sx q[0];
rz(-0.072364307) q[0];
x q[1];
rz(0.74283959) q[2];
sx q[2];
rz(-1.4577216) q[2];
sx q[2];
rz(-0.46087056) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0592812) q[1];
sx q[1];
rz(-2.0356405) q[1];
sx q[1];
rz(0.44537284) q[1];
x q[2];
rz(-0.71367587) q[3];
sx q[3];
rz(-1.6380966) q[3];
sx q[3];
rz(-2.0590559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7632873) q[2];
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
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0157938) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(0.75138599) q[0];
rz(1.3279351) q[1];
sx q[1];
rz(-1.8782047) q[1];
sx q[1];
rz(2.5352535) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060121814) q[0];
sx q[0];
rz(-1.6166286) q[0];
sx q[0];
rz(-1.5874552) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0688071) q[2];
sx q[2];
rz(-0.74337465) q[2];
sx q[2];
rz(0.547264) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1746754) q[1];
sx q[1];
rz(-1.7511909) q[1];
sx q[1];
rz(-1.0134407) q[1];
rz(-pi) q[2];
rz(1.8984406) q[3];
sx q[3];
rz(-0.22938211) q[3];
sx q[3];
rz(2.4250507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5027344) q[2];
sx q[2];
rz(-1.0417465) q[2];
sx q[2];
rz(1.139337) q[2];
rz(1.4849439) q[3];
sx q[3];
rz(-1.1805725) q[3];
sx q[3];
rz(3.0373354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6095603) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(-0.50810057) q[0];
rz(-1.5628901) q[1];
sx q[1];
rz(-1.0888313) q[1];
sx q[1];
rz(0.79024822) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3797487) q[0];
sx q[0];
rz(-1.2756691) q[0];
sx q[0];
rz(0.6167114) q[0];
rz(-pi) q[1];
rz(-0.054025606) q[2];
sx q[2];
rz(-0.21394357) q[2];
sx q[2];
rz(-0.33188785) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1192757) q[1];
sx q[1];
rz(-1.1806618) q[1];
sx q[1];
rz(-1.2783865) q[1];
rz(1.1931476) q[3];
sx q[3];
rz(-2.1834063) q[3];
sx q[3];
rz(0.8876422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.053085176) q[2];
sx q[2];
rz(-2.6997456) q[2];
sx q[2];
rz(-1.4132168) q[2];
rz(-1.4767856) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(3.0800381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4247894) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(2.0943663) q[0];
rz(2.5324902) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(-1.3887127) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6435218) q[0];
sx q[0];
rz(-2.5812491) q[0];
sx q[0];
rz(-1.6269006) q[0];
x q[1];
rz(1.9954761) q[2];
sx q[2];
rz(-1.0815902) q[2];
sx q[2];
rz(-2.3236772) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4365303) q[1];
sx q[1];
rz(-0.092518004) q[1];
sx q[1];
rz(-3.0180879) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8453313) q[3];
sx q[3];
rz(-0.98516609) q[3];
sx q[3];
rz(-1.1405917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1887112) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(-2.2593373) q[2];
rz(-1.4011718) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(1.2004948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27154487) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(-1.7154988) q[0];
rz(-3.0601314) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(0.55823278) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24960625) q[0];
sx q[0];
rz(-2.5476646) q[0];
sx q[0];
rz(-2.3114165) q[0];
rz(-2.0360557) q[2];
sx q[2];
rz(-0.12013398) q[2];
sx q[2];
rz(2.5533954) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8856814) q[1];
sx q[1];
rz(-2.3951888) q[1];
sx q[1];
rz(-1.4766775) q[1];
rz(-pi) q[2];
rz(2.8006323) q[3];
sx q[3];
rz(-0.51135671) q[3];
sx q[3];
rz(2.2380059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8490863) q[2];
sx q[2];
rz(-1.2693274) q[2];
sx q[2];
rz(-0.212184) q[2];
rz(2.9296181) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6367209) q[0];
sx q[0];
rz(-2.3175406) q[0];
sx q[0];
rz(-1.6037534) q[0];
rz(-2.3161855) q[1];
sx q[1];
rz(-2.4688265) q[1];
sx q[1];
rz(-2.6182981) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9742045) q[0];
sx q[0];
rz(-2.5295969) q[0];
sx q[0];
rz(-0.1083072) q[0];
rz(0.043838219) q[2];
sx q[2];
rz(-2.1423116) q[2];
sx q[2];
rz(2.8457355) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1112422) q[1];
sx q[1];
rz(-0.48644629) q[1];
sx q[1];
rz(0.72279795) q[1];
rz(-2.3862615) q[3];
sx q[3];
rz(-2.8108201) q[3];
sx q[3];
rz(-1.8473234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.837073) q[2];
sx q[2];
rz(-2.4261116) q[2];
sx q[2];
rz(0.26930299) q[2];
rz(-2.6473911) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(1.0860898) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8158648) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(1.6745463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(1.6124484) q[2];
sx q[2];
rz(-0.26298444) q[2];
sx q[2];
rz(2.0900805) q[2];
rz(0.79694637) q[3];
sx q[3];
rz(-1.4016101) q[3];
sx q[3];
rz(-2.3102643) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
