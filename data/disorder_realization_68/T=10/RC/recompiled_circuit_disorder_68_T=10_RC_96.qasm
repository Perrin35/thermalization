OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7473937) q[0];
sx q[0];
rz(-2.6497901) q[0];
sx q[0];
rz(2.9536182) q[0];
rz(-1.1176874) q[1];
sx q[1];
rz(-1.517065) q[1];
sx q[1];
rz(2.7741073) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7029019) q[0];
sx q[0];
rz(-0.54649788) q[0];
sx q[0];
rz(-1.370907) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8194524) q[2];
sx q[2];
rz(-2.6373632) q[2];
sx q[2];
rz(-2.3766975) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6403113) q[1];
sx q[1];
rz(-2.8456563) q[1];
sx q[1];
rz(-0.96247767) q[1];
rz(-pi) q[2];
rz(2.8005881) q[3];
sx q[3];
rz(-1.738428) q[3];
sx q[3];
rz(2.117702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1774896) q[2];
sx q[2];
rz(-2.6314645) q[2];
sx q[2];
rz(-2.5906079) q[2];
rz(1.8356813) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(1.3163542) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47857639) q[0];
sx q[0];
rz(-1.0000279) q[0];
sx q[0];
rz(0.4719032) q[0];
rz(0.42981237) q[1];
sx q[1];
rz(-1.8919573) q[1];
sx q[1];
rz(0.93634161) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82114007) q[0];
sx q[0];
rz(-1.3409412) q[0];
sx q[0];
rz(1.4810522) q[0];
rz(-pi) q[1];
rz(2.721644) q[2];
sx q[2];
rz(-2.311085) q[2];
sx q[2];
rz(-0.74479693) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97570005) q[1];
sx q[1];
rz(-0.70043889) q[1];
sx q[1];
rz(2.9794934) q[1];
x q[2];
rz(-1.068088) q[3];
sx q[3];
rz(-1.5503251) q[3];
sx q[3];
rz(2.3308144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3669746) q[2];
sx q[2];
rz(-0.32704157) q[2];
sx q[2];
rz(2.7152087) q[2];
rz(-1.2373699) q[3];
sx q[3];
rz(-0.62785134) q[3];
sx q[3];
rz(-3.1085076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8957829) q[0];
sx q[0];
rz(-1.824546) q[0];
sx q[0];
rz(2.202503) q[0];
rz(0.89871961) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(0.59392196) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9575189) q[0];
sx q[0];
rz(-0.44555095) q[0];
sx q[0];
rz(2.3222488) q[0];
x q[1];
rz(-0.8823231) q[2];
sx q[2];
rz(-1.5261298) q[2];
sx q[2];
rz(-0.67827144) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9582639) q[1];
sx q[1];
rz(-1.7839583) q[1];
sx q[1];
rz(2.0002685) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6731554) q[3];
sx q[3];
rz(-1.8803981) q[3];
sx q[3];
rz(-0.34164159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.64017355) q[2];
sx q[2];
rz(-2.4121425) q[2];
sx q[2];
rz(1.4397941) q[2];
rz(0.38763186) q[3];
sx q[3];
rz(-1.5165611) q[3];
sx q[3];
rz(0.38813996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3751635) q[0];
sx q[0];
rz(-1.5513993) q[0];
sx q[0];
rz(-0.50278062) q[0];
rz(0.76820961) q[1];
sx q[1];
rz(-0.50351024) q[1];
sx q[1];
rz(0.75685135) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3319791) q[0];
sx q[0];
rz(-0.50052128) q[0];
sx q[0];
rz(-0.74762263) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5966162) q[2];
sx q[2];
rz(-2.6757247) q[2];
sx q[2];
rz(0.70170882) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6824324) q[1];
sx q[1];
rz(-1.2130514) q[1];
sx q[1];
rz(2.1898502) q[1];
rz(2.4079011) q[3];
sx q[3];
rz(-1.1847704) q[3];
sx q[3];
rz(1.8611849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42671529) q[2];
sx q[2];
rz(-1.2299512) q[2];
sx q[2];
rz(1.4871917) q[2];
rz(-2.5590844) q[3];
sx q[3];
rz(-1.094386) q[3];
sx q[3];
rz(-2.5845161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29397598) q[0];
sx q[0];
rz(-2.0563545) q[0];
sx q[0];
rz(-2.3838682) q[0];
rz(-1.853653) q[1];
sx q[1];
rz(-0.92823354) q[1];
sx q[1];
rz(2.0910738) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61705631) q[0];
sx q[0];
rz(-2.7334088) q[0];
sx q[0];
rz(0.88390669) q[0];
rz(-pi) q[1];
rz(-0.96111416) q[2];
sx q[2];
rz(-2.4498307) q[2];
sx q[2];
rz(1.320653) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1615636) q[1];
sx q[1];
rz(-2.1254351) q[1];
sx q[1];
rz(0.6273004) q[1];
x q[2];
rz(-1.0195144) q[3];
sx q[3];
rz(-0.84421221) q[3];
sx q[3];
rz(-0.36908484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.22333764) q[2];
sx q[2];
rz(-2.788322) q[2];
sx q[2];
rz(-2.5081432) q[2];
rz(1.9472306) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(-0.71715322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69960064) q[0];
sx q[0];
rz(-2.6514335) q[0];
sx q[0];
rz(2.8884086) q[0];
rz(1.5340012) q[1];
sx q[1];
rz(-1.7065159) q[1];
sx q[1];
rz(-1.6794499) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40186858) q[0];
sx q[0];
rz(-1.7812294) q[0];
sx q[0];
rz(-1.6519288) q[0];
x q[1];
rz(-0.24121933) q[2];
sx q[2];
rz(-2.2058645) q[2];
sx q[2];
rz(1.5347753) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6219382) q[1];
sx q[1];
rz(-2.7192273) q[1];
sx q[1];
rz(-2.7154891) q[1];
rz(-pi) q[2];
rz(2.3092689) q[3];
sx q[3];
rz(-0.451085) q[3];
sx q[3];
rz(-3.0576599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9399461) q[2];
sx q[2];
rz(-0.74735171) q[2];
sx q[2];
rz(0.80491006) q[2];
rz(-1.1770052) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(2.0578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.8686304) q[0];
sx q[0];
rz(-1.0694163) q[0];
sx q[0];
rz(0.92765635) q[0];
rz(-2.1169128) q[1];
sx q[1];
rz(-1.6352446) q[1];
sx q[1];
rz(2.129508) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0439118) q[0];
sx q[0];
rz(-1.0543531) q[0];
sx q[0];
rz(-0.54215706) q[0];
rz(-1.7094678) q[2];
sx q[2];
rz(-1.1853293) q[2];
sx q[2];
rz(-3.1388381) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.583657) q[1];
sx q[1];
rz(-2.5968938) q[1];
sx q[1];
rz(-3.0779561) q[1];
x q[2];
rz(-0.55862553) q[3];
sx q[3];
rz(-1.6405676) q[3];
sx q[3];
rz(-0.39486265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6107789) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(-0.42993316) q[2];
rz(-2.1271465) q[3];
sx q[3];
rz(-2.7323664) q[3];
sx q[3];
rz(-2.608192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5291418) q[0];
sx q[0];
rz(-2.2021459) q[0];
sx q[0];
rz(2.9274143) q[0];
rz(-1.051349) q[1];
sx q[1];
rz(-0.21251692) q[1];
sx q[1];
rz(2.8578551) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34127125) q[0];
sx q[0];
rz(-0.30297908) q[0];
sx q[0];
rz(3.0269701) q[0];
x q[1];
rz(-1.0546513) q[2];
sx q[2];
rz(-0.37545855) q[2];
sx q[2];
rz(-1.5309389) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.46854308) q[1];
sx q[1];
rz(-1.0769516) q[1];
sx q[1];
rz(1.7935351) q[1];
rz(-pi) q[2];
rz(0.78332087) q[3];
sx q[3];
rz(-2.3187227) q[3];
sx q[3];
rz(-2.2058723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.45067898) q[2];
sx q[2];
rz(-2.6288855) q[2];
sx q[2];
rz(1.4452176) q[2];
rz(1.5444267) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(-2.8022695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.504869) q[0];
sx q[0];
rz(-0.050014194) q[0];
sx q[0];
rz(-0.069256393) q[0];
rz(1.4878558) q[1];
sx q[1];
rz(-1.2830877) q[1];
sx q[1];
rz(1.5725296) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3587787) q[0];
sx q[0];
rz(-2.1920715) q[0];
sx q[0];
rz(-0.57399477) q[0];
rz(-pi) q[1];
rz(1.6530767) q[2];
sx q[2];
rz(-1.8956208) q[2];
sx q[2];
rz(-0.61330739) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.931817) q[1];
sx q[1];
rz(-1.4869542) q[1];
sx q[1];
rz(0.2340338) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2980372) q[3];
sx q[3];
rz(-2.4747304) q[3];
sx q[3];
rz(0.71344261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9562324) q[2];
sx q[2];
rz(-2.9113443) q[2];
sx q[2];
rz(-3.0017079) q[2];
rz(-2.774003) q[3];
sx q[3];
rz(-1.1871754) q[3];
sx q[3];
rz(0.99115133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763828) q[0];
sx q[0];
rz(-0.39127025) q[0];
sx q[0];
rz(2.4998253) q[0];
rz(1.2311252) q[1];
sx q[1];
rz(-1.9893913) q[1];
sx q[1];
rz(0.26783255) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56775996) q[0];
sx q[0];
rz(-2.0127675) q[0];
sx q[0];
rz(2.3964336) q[0];
rz(-pi) q[1];
rz(-1.3875302) q[2];
sx q[2];
rz(-1.4472618) q[2];
sx q[2];
rz(-1.087041) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6707582) q[1];
sx q[1];
rz(-1.8948312) q[1];
sx q[1];
rz(1.8761937) q[1];
x q[2];
rz(-1.1166499) q[3];
sx q[3];
rz(-0.59026679) q[3];
sx q[3];
rz(-0.86436194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3315167) q[2];
sx q[2];
rz(-1.5248359) q[2];
sx q[2];
rz(-1.4769185) q[2];
rz(-2.8752575) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(0.5464856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01263604) q[0];
sx q[0];
rz(-2.2205882) q[0];
sx q[0];
rz(0.90482774) q[0];
rz(-2.3616882) q[1];
sx q[1];
rz(-2.654568) q[1];
sx q[1];
rz(1.754896) q[1];
rz(0.65200381) q[2];
sx q[2];
rz(-0.82521385) q[2];
sx q[2];
rz(3.0575976) q[2];
rz(-1.3569309) q[3];
sx q[3];
rz(-2.2371117) q[3];
sx q[3];
rz(-0.37123751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
