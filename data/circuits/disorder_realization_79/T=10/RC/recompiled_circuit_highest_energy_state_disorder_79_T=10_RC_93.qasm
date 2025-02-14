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
rz(-2.7231554) q[0];
sx q[0];
rz(-2.1783481) q[0];
sx q[0];
rz(2.9377687) q[0];
rz(1.0526429) q[1];
sx q[1];
rz(3.9349603) q[1];
sx q[1];
rz(10.959672) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2685854) q[0];
sx q[0];
rz(-2.632336) q[0];
sx q[0];
rz(-1.0727706) q[0];
rz(-2.7947172) q[2];
sx q[2];
rz(-0.87793186) q[2];
sx q[2];
rz(1.3077867) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.71434778) q[1];
sx q[1];
rz(-0.74810076) q[1];
sx q[1];
rz(2.1022878) q[1];
rz(-0.20333692) q[3];
sx q[3];
rz(-1.4383738) q[3];
sx q[3];
rz(1.9258969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7242929) q[2];
sx q[2];
rz(-2.5827926) q[2];
sx q[2];
rz(-0.99754769) q[2];
rz(-0.7558465) q[3];
sx q[3];
rz(-1.7852802) q[3];
sx q[3];
rz(-2.1006179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.7928829) q[0];
sx q[0];
rz(-3.0443158) q[0];
sx q[0];
rz(1.2874228) q[0];
rz(-0.48311326) q[1];
sx q[1];
rz(-1.0419507) q[1];
sx q[1];
rz(1.2189254) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.481015) q[0];
sx q[0];
rz(-1.7661375) q[0];
sx q[0];
rz(-0.16150772) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4697815) q[2];
sx q[2];
rz(-1.3992568) q[2];
sx q[2];
rz(-1.1900657) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9049553) q[1];
sx q[1];
rz(-1.4262916) q[1];
sx q[1];
rz(1.2895686) q[1];
x q[2];
rz(1.990149) q[3];
sx q[3];
rz(-2.4046728) q[3];
sx q[3];
rz(1.9896955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7634742) q[2];
sx q[2];
rz(-3.0749622) q[2];
sx q[2];
rz(0.62180579) q[2];
rz(2.5032737) q[3];
sx q[3];
rz(-0.74898762) q[3];
sx q[3];
rz(-0.84426713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3062375) q[0];
sx q[0];
rz(-0.64202809) q[0];
sx q[0];
rz(-2.9252692) q[0];
rz(0.59858876) q[1];
sx q[1];
rz(-1.8363771) q[1];
sx q[1];
rz(-2.6187706) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3591839) q[0];
sx q[0];
rz(-2.6893927) q[0];
sx q[0];
rz(0.81069273) q[0];
rz(-0.39117809) q[2];
sx q[2];
rz(-1.7157946) q[2];
sx q[2];
rz(-0.36996335) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.93980689) q[1];
sx q[1];
rz(-2.5720189) q[1];
sx q[1];
rz(-2.7199183) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9013635) q[3];
sx q[3];
rz(-1.2798314) q[3];
sx q[3];
rz(1.3823079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9598976) q[2];
sx q[2];
rz(-1.8671702) q[2];
sx q[2];
rz(1.5175021) q[2];
rz(0.46487871) q[3];
sx q[3];
rz(-0.93022323) q[3];
sx q[3];
rz(-1.5215065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14438039) q[0];
sx q[0];
rz(-2.7535487) q[0];
sx q[0];
rz(-1.5390747) q[0];
rz(-1.0333565) q[1];
sx q[1];
rz(-0.76834232) q[1];
sx q[1];
rz(1.8446406) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42596969) q[0];
sx q[0];
rz(-2.2821626) q[0];
sx q[0];
rz(2.864896) q[0];
rz(0.06242604) q[2];
sx q[2];
rz(-0.87597825) q[2];
sx q[2];
rz(-2.7165987) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0721283) q[1];
sx q[1];
rz(-1.1934682) q[1];
sx q[1];
rz(1.1262141) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55418535) q[3];
sx q[3];
rz(-1.300214) q[3];
sx q[3];
rz(0.24662805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3492655) q[2];
sx q[2];
rz(-2.5040864) q[2];
sx q[2];
rz(-2.0959334) q[2];
rz(-0.21444923) q[3];
sx q[3];
rz(-2.4172473) q[3];
sx q[3];
rz(-1.9946056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3359208) q[0];
sx q[0];
rz(-1.2821953) q[0];
sx q[0];
rz(-0.51934284) q[0];
rz(-2.1995811) q[1];
sx q[1];
rz(-2.8603034) q[1];
sx q[1];
rz(1.5325783) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62270861) q[0];
sx q[0];
rz(-1.2075338) q[0];
sx q[0];
rz(-2.4469482) q[0];
rz(-1.4283044) q[2];
sx q[2];
rz(-1.7497369) q[2];
sx q[2];
rz(-1.6662836) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6841507) q[1];
sx q[1];
rz(-1.5655715) q[1];
sx q[1];
rz(-2.0713776) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4761837) q[3];
sx q[3];
rz(-1.4615371) q[3];
sx q[3];
rz(-0.13302468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4970826) q[2];
sx q[2];
rz(-2.3641059) q[2];
sx q[2];
rz(-0.28437781) q[2];
rz(0.55822462) q[3];
sx q[3];
rz(-1.3642718) q[3];
sx q[3];
rz(-2.8936581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.072642) q[0];
sx q[0];
rz(-2.729029) q[0];
sx q[0];
rz(0.64055842) q[0];
rz(-1.6259646) q[1];
sx q[1];
rz(-1.1311572) q[1];
sx q[1];
rz(0.21387771) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.36422) q[0];
sx q[0];
rz(-1.3427444) q[0];
sx q[0];
rz(1.8814726) q[0];
rz(-0.74486129) q[2];
sx q[2];
rz(-1.8554885) q[2];
sx q[2];
rz(-2.8887951) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2474413) q[1];
sx q[1];
rz(-0.94235984) q[1];
sx q[1];
rz(3.0559866) q[1];
x q[2];
rz(1.6951896) q[3];
sx q[3];
rz(-0.52342192) q[3];
sx q[3];
rz(2.058213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3512257) q[2];
sx q[2];
rz(-2.6595317) q[2];
sx q[2];
rz(-2.9638929) q[2];
rz(-1.7443582) q[3];
sx q[3];
rz(-1.9652365) q[3];
sx q[3];
rz(2.7241404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099667065) q[0];
sx q[0];
rz(-0.61138994) q[0];
sx q[0];
rz(-2.0360816) q[0];
rz(2.8583177) q[1];
sx q[1];
rz(-2.6073644) q[1];
sx q[1];
rz(-1.4615321) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7992226) q[0];
sx q[0];
rz(-2.2307751) q[0];
sx q[0];
rz(-2.2837385) q[0];
rz(-2.6335476) q[2];
sx q[2];
rz(-0.99544243) q[2];
sx q[2];
rz(-2.3066556) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4883621) q[1];
sx q[1];
rz(-1.686859) q[1];
sx q[1];
rz(-0.56452063) q[1];
rz(-pi) q[2];
rz(1.8959778) q[3];
sx q[3];
rz(-2.1228613) q[3];
sx q[3];
rz(-2.0887041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0561515) q[2];
sx q[2];
rz(-1.5957007) q[2];
sx q[2];
rz(2.936787) q[2];
rz(1.0497931) q[3];
sx q[3];
rz(-0.27725163) q[3];
sx q[3];
rz(-2.7753196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3664704) q[0];
sx q[0];
rz(-2.6752052) q[0];
sx q[0];
rz(-3.108016) q[0];
rz(-1.6665943) q[1];
sx q[1];
rz(-0.74445236) q[1];
sx q[1];
rz(1.9974476) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1468723) q[0];
sx q[0];
rz(-2.0573186) q[0];
sx q[0];
rz(-1.542132) q[0];
rz(-1.894925) q[2];
sx q[2];
rz(-1.4767382) q[2];
sx q[2];
rz(-0.45378527) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5976286) q[1];
sx q[1];
rz(-1.9457327) q[1];
sx q[1];
rz(0.72245325) q[1];
rz(-pi) q[2];
x q[2];
rz(2.860922) q[3];
sx q[3];
rz(-1.2378527) q[3];
sx q[3];
rz(1.5047634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7555776) q[2];
sx q[2];
rz(-0.80499804) q[2];
sx q[2];
rz(-1.8899274) q[2];
rz(1.3517316) q[3];
sx q[3];
rz(-0.91910619) q[3];
sx q[3];
rz(-2.1613817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1397454) q[0];
sx q[0];
rz(-2.503105) q[0];
sx q[0];
rz(1.8852604) q[0];
rz(-0.095257692) q[1];
sx q[1];
rz(-0.88905159) q[1];
sx q[1];
rz(-2.6108066) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4149218) q[0];
sx q[0];
rz(-1.7309395) q[0];
sx q[0];
rz(-1.3199575) q[0];
rz(2.8235648) q[2];
sx q[2];
rz(-1.6075168) q[2];
sx q[2];
rz(-1.5098234) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.8618882) q[1];
sx q[1];
rz(-1.8100396) q[1];
sx q[1];
rz(-1.4390535) q[1];
x q[2];
rz(1.68566) q[3];
sx q[3];
rz(-1.126976) q[3];
sx q[3];
rz(2.165909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.329616) q[2];
sx q[2];
rz(-1.5584385) q[2];
sx q[2];
rz(-2.6007593) q[2];
rz(2.3234308) q[3];
sx q[3];
rz(-1.6988138) q[3];
sx q[3];
rz(2.5087859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.624991) q[0];
sx q[0];
rz(-1.7778968) q[0];
sx q[0];
rz(2.7428108) q[0];
rz(-1.3840236) q[1];
sx q[1];
rz(-0.75474352) q[1];
sx q[1];
rz(-3.0922281) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20959768) q[0];
sx q[0];
rz(-1.6315797) q[0];
sx q[0];
rz(1.6796735) q[0];
x q[1];
rz(0.10714679) q[2];
sx q[2];
rz(-1.309589) q[2];
sx q[2];
rz(3.0334453) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70061848) q[1];
sx q[1];
rz(-1.5262394) q[1];
sx q[1];
rz(0.99127165) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3672164) q[3];
sx q[3];
rz(-0.94881159) q[3];
sx q[3];
rz(-0.79727117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0177239) q[2];
sx q[2];
rz(-1.3088635) q[2];
sx q[2];
rz(-1.6309942) q[2];
rz(-2.2137568) q[3];
sx q[3];
rz(-1.8001013) q[3];
sx q[3];
rz(1.9063037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8265726) q[0];
sx q[0];
rz(-2.6980504) q[0];
sx q[0];
rz(2.0499688) q[0];
rz(-2.0768968) q[1];
sx q[1];
rz(-2.6152492) q[1];
sx q[1];
rz(-1.1660887) q[1];
rz(1.0749558) q[2];
sx q[2];
rz(-1.1100162) q[2];
sx q[2];
rz(0.77592862) q[2];
rz(-0.31872411) q[3];
sx q[3];
rz(-2.095185) q[3];
sx q[3];
rz(-0.17879055) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
