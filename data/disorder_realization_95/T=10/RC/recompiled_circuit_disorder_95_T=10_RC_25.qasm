OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.99958217) q[0];
sx q[0];
rz(5.2810623) q[0];
sx q[0];
rz(5.3856344) q[0];
rz(-0.23437962) q[1];
sx q[1];
rz(-0.27581629) q[1];
sx q[1];
rz(2.0770567) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9165186) q[0];
sx q[0];
rz(-1.5009891) q[0];
sx q[0];
rz(-2.201447) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2034982) q[2];
sx q[2];
rz(-1.2812867) q[2];
sx q[2];
rz(-3.0245568) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0734288) q[1];
sx q[1];
rz(-2.4792719) q[1];
sx q[1];
rz(-2.8627002) q[1];
x q[2];
rz(0.75818054) q[3];
sx q[3];
rz(-1.7054134) q[3];
sx q[3];
rz(1.7077703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.477318) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(2.409639) q[2];
rz(2.1814363) q[3];
sx q[3];
rz(-0.82257661) q[3];
sx q[3];
rz(1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17523781) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(-0.91645855) q[0];
rz(0.48049277) q[1];
sx q[1];
rz(-2.5669211) q[1];
sx q[1];
rz(2.2629471) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4247596) q[0];
sx q[0];
rz(-0.75100198) q[0];
sx q[0];
rz(-0.99320937) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8354561) q[2];
sx q[2];
rz(-2.172643) q[2];
sx q[2];
rz(-2.2167609) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.52777427) q[1];
sx q[1];
rz(-1.4973745) q[1];
sx q[1];
rz(0.41927494) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.068275498) q[3];
sx q[3];
rz(-1.8797415) q[3];
sx q[3];
rz(0.49648778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2800704) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(2.4943165) q[2];
rz(-2.9679126) q[3];
sx q[3];
rz(-2.0208385) q[3];
sx q[3];
rz(0.15163264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.52755255) q[0];
sx q[0];
rz(-2.5019167) q[0];
sx q[0];
rz(-2.8955984) q[0];
rz(1.4100769) q[1];
sx q[1];
rz(-1.9672111) q[1];
sx q[1];
rz(2.0203967) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4778053) q[0];
sx q[0];
rz(-0.87513798) q[0];
sx q[0];
rz(1.0870766) q[0];
x q[1];
rz(-1.4175225) q[2];
sx q[2];
rz(-0.3970662) q[2];
sx q[2];
rz(3.0561471) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.4157297) q[1];
sx q[1];
rz(-1.2408857) q[1];
sx q[1];
rz(1.2582474) q[1];
rz(2.9259053) q[3];
sx q[3];
rz(-2.6378999) q[3];
sx q[3];
rz(1.4195201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7401509) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(0.95139727) q[2];
rz(-0.65008632) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(-2.8459809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6275416) q[0];
sx q[0];
rz(-2.5294332) q[0];
sx q[0];
rz(2.3535368) q[0];
rz(1.0568985) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(-0.11638164) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48556337) q[0];
sx q[0];
rz(-2.2040327) q[0];
sx q[0];
rz(-1.5390736) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1999947) q[2];
sx q[2];
rz(-2.5430352) q[2];
sx q[2];
rz(1.7184005) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1981922) q[1];
sx q[1];
rz(-2.0640089) q[1];
sx q[1];
rz(-2.0342159) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0451123) q[3];
sx q[3];
rz(-2.6498142) q[3];
sx q[3];
rz(-1.5887367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5300166) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(0.25203618) q[2];
rz(-2.7633372) q[3];
sx q[3];
rz(-0.16246048) q[3];
sx q[3];
rz(2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56548059) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(-2.6089923) q[0];
rz(-1.700092) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(1.1486357) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.900433) q[0];
sx q[0];
rz(-0.94928375) q[0];
sx q[0];
rz(-0.09089367) q[0];
x q[1];
rz(2.1158475) q[2];
sx q[2];
rz(-2.3902241) q[2];
sx q[2];
rz(2.2124706) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.83541218) q[1];
sx q[1];
rz(-0.40459834) q[1];
sx q[1];
rz(-0.49680357) q[1];
rz(-1.1762189) q[3];
sx q[3];
rz(-0.88486457) q[3];
sx q[3];
rz(-1.8358313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1308412) q[2];
sx q[2];
rz(-0.66117078) q[2];
sx q[2];
rz(-1.032069) q[2];
rz(0.71470913) q[3];
sx q[3];
rz(-1.8604449) q[3];
sx q[3];
rz(-0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.5500568) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(-0.39598879) q[0];
rz(1.6962601) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(-2.9352303) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43564046) q[0];
sx q[0];
rz(-0.75169509) q[0];
sx q[0];
rz(-0.12775001) q[0];
rz(-0.80870734) q[2];
sx q[2];
rz(-2.2924097) q[2];
sx q[2];
rz(-1.8546113) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5620835) q[1];
sx q[1];
rz(-1.675203) q[1];
sx q[1];
rz(3.0444006) q[1];
x q[2];
rz(-1.0422802) q[3];
sx q[3];
rz(-2.7665666) q[3];
sx q[3];
rz(0.11881766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6017194) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(-2.8708141) q[2];
rz(2.3932636) q[3];
sx q[3];
rz(-1.7691408) q[3];
sx q[3];
rz(1.07871) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2449743) q[0];
sx q[0];
rz(-2.0386319) q[0];
sx q[0];
rz(-0.57624972) q[0];
rz(2.9684864) q[1];
sx q[1];
rz(-2.3997967) q[1];
sx q[1];
rz(1.9304088) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9528708) q[0];
sx q[0];
rz(-1.0591918) q[0];
sx q[0];
rz(-2.1923724) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7089825) q[2];
sx q[2];
rz(-1.0850564) q[2];
sx q[2];
rz(-1.1952343) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.83963517) q[1];
sx q[1];
rz(-0.17050276) q[1];
sx q[1];
rz(-1.9819928) q[1];
rz(-pi) q[2];
rz(2.1738449) q[3];
sx q[3];
rz(-1.839404) q[3];
sx q[3];
rz(-0.3097765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8048191) q[2];
sx q[2];
rz(-2.4833198) q[2];
sx q[2];
rz(-0.024519196) q[2];
rz(-2.426614) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(-1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1361168) q[0];
sx q[0];
rz(-1.5908717) q[0];
sx q[0];
rz(1.0205644) q[0];
rz(0.15469805) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(1.9205836) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8120136) q[0];
sx q[0];
rz(-2.3276969) q[0];
sx q[0];
rz(1.194186) q[0];
rz(-pi) q[1];
rz(1.4079354) q[2];
sx q[2];
rz(-1.3961627) q[2];
sx q[2];
rz(-2.3101431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.11215969) q[1];
sx q[1];
rz(-2.9661848) q[1];
sx q[1];
rz(2.4584241) q[1];
rz(-pi) q[2];
rz(0.82151316) q[3];
sx q[3];
rz(-1.8213846) q[3];
sx q[3];
rz(-2.487395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.28875479) q[2];
sx q[2];
rz(-1.6483665) q[2];
sx q[2];
rz(2.4364046) q[2];
rz(0.60338902) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(-1.6220629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.110638) q[0];
sx q[0];
rz(-1.5663261) q[0];
sx q[0];
rz(-0.13701339) q[0];
rz(-2.5367472) q[1];
sx q[1];
rz(-0.73692656) q[1];
sx q[1];
rz(3.0158214) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.754697) q[0];
sx q[0];
rz(-1.7594975) q[0];
sx q[0];
rz(-0.17908355) q[0];
rz(-pi) q[1];
rz(0.33340402) q[2];
sx q[2];
rz(-1.3840904) q[2];
sx q[2];
rz(-0.29078996) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3198493) q[1];
sx q[1];
rz(-0.84801596) q[1];
sx q[1];
rz(-2.7601932) q[1];
rz(1.3689234) q[3];
sx q[3];
rz(-1.2807506) q[3];
sx q[3];
rz(1.2253075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.003309) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(2.4323145) q[2];
rz(-2.5214031) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(2.9848849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.9897292) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(-1.9412769) q[0];
rz(-0.26750803) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(-1.1402003) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32894293) q[0];
sx q[0];
rz(-0.2812627) q[0];
sx q[0];
rz(1.9360696) q[0];
rz(2.0613725) q[2];
sx q[2];
rz(-0.51418257) q[2];
sx q[2];
rz(2.0928004) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0298437) q[1];
sx q[1];
rz(-0.20856253) q[1];
sx q[1];
rz(-0.44058056) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8885918) q[3];
sx q[3];
rz(-0.80298775) q[3];
sx q[3];
rz(0.12249891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.37529477) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(-0.89861384) q[2];
rz(-0.0071772655) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(-1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89467775) q[0];
sx q[0];
rz(-1.4151731) q[0];
sx q[0];
rz(-1.7807501) q[0];
rz(0.5207516) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(-0.76277914) q[2];
sx q[2];
rz(-0.63763466) q[2];
sx q[2];
rz(0.55479738) q[2];
rz(-0.58581523) q[3];
sx q[3];
rz(-1.9184434) q[3];
sx q[3];
rz(-0.8003269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
