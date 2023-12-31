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
rz(0.51529348) q[0];
rz(1.3442858) q[1];
sx q[1];
rz(-2.9872515) q[1];
sx q[1];
rz(0.57758346) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2429457) q[0];
sx q[0];
rz(-2.4624914) q[0];
sx q[0];
rz(-0.26429096) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53519997) q[2];
sx q[2];
rz(-2.0173965) q[2];
sx q[2];
rz(-1.5314147) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9211728) q[1];
sx q[1];
rz(-2.6487659) q[1];
sx q[1];
rz(-0.9763896) q[1];
rz(-pi) q[2];
rz(0.9382117) q[3];
sx q[3];
rz(-1.0217561) q[3];
sx q[3];
rz(0.051018056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43705964) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(-2.4543767) q[2];
rz(2.1263188) q[3];
sx q[3];
rz(-1.7679368) q[3];
sx q[3];
rz(-3.0190873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9706443) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(1.2600391) q[0];
rz(-2.1353703) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(0.84567436) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0802404) q[0];
sx q[0];
rz(-2.4363359) q[0];
sx q[0];
rz(-2.4387226) q[0];
rz(-pi) q[1];
rz(-3.0078366) q[2];
sx q[2];
rz(-1.4417366) q[2];
sx q[2];
rz(-1.2226906) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0963124) q[1];
sx q[1];
rz(-1.7573866) q[1];
sx q[1];
rz(-2.6394096) q[1];
rz(-pi) q[2];
rz(-0.50176974) q[3];
sx q[3];
rz(-2.1309149) q[3];
sx q[3];
rz(-3.113941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.78850293) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(-2.6611924) q[2];
rz(1.3530312) q[3];
sx q[3];
rz(-1.055911) q[3];
sx q[3];
rz(1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1903494) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(-2.3685266) q[0];
rz(0.13126016) q[1];
sx q[1];
rz(-1.8570329) q[1];
sx q[1];
rz(-2.0551596) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4341136) q[0];
sx q[0];
rz(-1.1093603) q[0];
sx q[0];
rz(-3.1397318) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92347446) q[2];
sx q[2];
rz(-1.662865) q[2];
sx q[2];
rz(2.0331241) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74124386) q[1];
sx q[1];
rz(-1.8596008) q[1];
sx q[1];
rz(1.5762394) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9757189) q[3];
sx q[3];
rz(-1.4110663) q[3];
sx q[3];
rz(-0.83042849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0125668) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(1.770299) q[2];
rz(-0.38315547) q[3];
sx q[3];
rz(-1.8846735) q[3];
sx q[3];
rz(-2.3390521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4521769) q[0];
sx q[0];
rz(-1.2503662) q[0];
sx q[0];
rz(-3.0932328) q[0];
rz(0.16391779) q[1];
sx q[1];
rz(-0.36968958) q[1];
sx q[1];
rz(1.4455459) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3460798) q[0];
sx q[0];
rz(-2.4582986) q[0];
sx q[0];
rz(2.8542095) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7518026) q[2];
sx q[2];
rz(-0.23332694) q[2];
sx q[2];
rz(-2.0248272) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2913937) q[1];
sx q[1];
rz(-1.7899917) q[1];
sx q[1];
rz(1.5047969) q[1];
x q[2];
rz(-1.0936071) q[3];
sx q[3];
rz(-0.573199) q[3];
sx q[3];
rz(-1.5968061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.066102862) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(2.1172822) q[2];
rz(-1.6131489) q[3];
sx q[3];
rz(-1.5214835) q[3];
sx q[3];
rz(2.9083692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7016474) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(1.8540927) q[0];
rz(-2.8225186) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(-0.85420001) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4860977) q[0];
sx q[0];
rz(-0.71394074) q[0];
sx q[0];
rz(1.5582725) q[0];
x q[1];
rz(2.6536077) q[2];
sx q[2];
rz(-0.93163604) q[2];
sx q[2];
rz(2.0069063) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8728719) q[1];
sx q[1];
rz(-0.99209058) q[1];
sx q[1];
rz(2.242356) q[1];
rz(-pi) q[2];
rz(-0.79980005) q[3];
sx q[3];
rz(-2.4197257) q[3];
sx q[3];
rz(2.0197899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.95191082) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(2.5642776) q[2];
rz(0.50950766) q[3];
sx q[3];
rz(-2.7039492) q[3];
sx q[3];
rz(-0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1381056) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(-3.0694718) q[0];
rz(2.0347319) q[1];
sx q[1];
rz(-0.51263428) q[1];
sx q[1];
rz(-3.0153826) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.431625) q[0];
sx q[0];
rz(-0.21490782) q[0];
sx q[0];
rz(-2.9931195) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1926786) q[2];
sx q[2];
rz(-0.83231229) q[2];
sx q[2];
rz(-0.2371012) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.63657657) q[1];
sx q[1];
rz(-2.2762183) q[1];
sx q[1];
rz(-0.24800639) q[1];
rz(-pi) q[2];
rz(-2.0893655) q[3];
sx q[3];
rz(-0.37399451) q[3];
sx q[3];
rz(2.3911589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2386027) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(0.77511707) q[2];
rz(-2.3136247) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(2.0194139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59259748) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(3.0431842) q[0];
rz(1.9495643) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(-0.55955204) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0347621) q[0];
sx q[0];
rz(-2.1666399) q[0];
sx q[0];
rz(-1.725127) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6227116) q[2];
sx q[2];
rz(-1.8318818) q[2];
sx q[2];
rz(-1.7920997) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65161381) q[1];
sx q[1];
rz(-1.4414806) q[1];
sx q[1];
rz(3.053385) q[1];
rz(-pi) q[2];
rz(-1.7667174) q[3];
sx q[3];
rz(-2.19176) q[3];
sx q[3];
rz(-1.8718815) q[3];
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
rz(-0.5665468) q[3];
sx q[3];
rz(-2.6930801) q[3];
sx q[3];
rz(-0.47376537) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.375181) q[0];
sx q[0];
rz(-1.8090929) q[0];
sx q[0];
rz(0.73076105) q[0];
rz(-0.14239755) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(-2.2699845) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7155834) q[0];
sx q[0];
rz(-1.6462438) q[0];
sx q[0];
rz(1.5619318) q[0];
rz(-pi) q[1];
rz(2.7996054) q[2];
sx q[2];
rz(-2.7139398) q[2];
sx q[2];
rz(0.74117408) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3155047) q[1];
sx q[1];
rz(-1.9313889) q[1];
sx q[1];
rz(-2.4396067) q[1];
rz(-0.96846795) q[3];
sx q[3];
rz(-1.6864711) q[3];
sx q[3];
rz(-2.065421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4153851) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(1.139572) q[2];
rz(1.6428044) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(2.22877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20294872) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(-1.9198445) q[0];
rz(-0.16601673) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(-1.5244012) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071951889) q[0];
sx q[0];
rz(-1.4884293) q[0];
sx q[0];
rz(-0.029043341) q[0];
rz(-pi) q[1];
rz(-1.7622856) q[2];
sx q[2];
rz(-2.0404437) q[2];
sx q[2];
rz(-0.42275235) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0601378) q[1];
sx q[1];
rz(-0.36768915) q[1];
sx q[1];
rz(1.1170438) q[1];
rz(-pi) q[2];
rz(-1.2576305) q[3];
sx q[3];
rz(-2.1861665) q[3];
sx q[3];
rz(0.19389158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4979424) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(-1.3051916) q[0];
rz(-0.38048831) q[1];
sx q[1];
rz(-2.0996129) q[1];
sx q[1];
rz(-0.25340733) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6127374) q[0];
sx q[0];
rz(-2.7326072) q[0];
sx q[0];
rz(-1.3091875) q[0];
rz(1.2316522) q[2];
sx q[2];
rz(-1.7288496) q[2];
sx q[2];
rz(0.16976419) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39292654) q[1];
sx q[1];
rz(-1.8112507) q[1];
sx q[1];
rz(-2.7886831) q[1];
rz(-0.6110544) q[3];
sx q[3];
rz(-1.9285893) q[3];
sx q[3];
rz(2.4886481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5499251) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(-2.6386476) q[2];
rz(-0.89899603) q[3];
sx q[3];
rz(-1.2939724) q[3];
sx q[3];
rz(-1.1635273) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1702561) q[0];
sx q[0];
rz(-1.6032871) q[0];
sx q[0];
rz(0.26300318) q[0];
rz(-0.7111711) q[1];
sx q[1];
rz(-1.0881337) q[1];
sx q[1];
rz(1.7137391) q[1];
rz(2.860643) q[2];
sx q[2];
rz(-1.3215669) q[2];
sx q[2];
rz(2.486698) q[2];
rz(2.4651299) q[3];
sx q[3];
rz(-2.2623895) q[3];
sx q[3];
rz(-1.1415979) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
