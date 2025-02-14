OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1333753) q[0];
sx q[0];
rz(-1.8741338) q[0];
sx q[0];
rz(0.01292364) q[0];
rz(-2.456993) q[1];
sx q[1];
rz(-0.79889387) q[1];
sx q[1];
rz(2.0838783) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3290214) q[0];
sx q[0];
rz(-2.2264535) q[0];
sx q[0];
rz(1.283487) q[0];
rz(-2.3637412) q[2];
sx q[2];
rz(-1.9490644) q[2];
sx q[2];
rz(-0.08535484) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5200708) q[1];
sx q[1];
rz(-2.3316262) q[1];
sx q[1];
rz(2.6317627) q[1];
rz(-pi) q[2];
rz(3.0989981) q[3];
sx q[3];
rz(-1.179368) q[3];
sx q[3];
rz(-0.0022050641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5554123) q[2];
sx q[2];
rz(-1.5427898) q[2];
sx q[2];
rz(3.0482698) q[2];
rz(0.020545067) q[3];
sx q[3];
rz(-2.8829657) q[3];
sx q[3];
rz(-1.3625905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0697295) q[0];
sx q[0];
rz(-1.7427895) q[0];
sx q[0];
rz(-2.3139957) q[0];
rz(-2.1780275) q[1];
sx q[1];
rz(-1.6160485) q[1];
sx q[1];
rz(-0.73659426) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73294357) q[0];
sx q[0];
rz(-1.8445065) q[0];
sx q[0];
rz(-3.0628171) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57666333) q[2];
sx q[2];
rz(-2.0915589) q[2];
sx q[2];
rz(1.3530089) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0328278) q[1];
sx q[1];
rz(-1.2695644) q[1];
sx q[1];
rz(2.9297622) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16421825) q[3];
sx q[3];
rz(-0.13623691) q[3];
sx q[3];
rz(3.1395903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5654512) q[2];
sx q[2];
rz(-2.5665923) q[2];
sx q[2];
rz(2.3973993) q[2];
rz(-0.60892504) q[3];
sx q[3];
rz(-0.78151339) q[3];
sx q[3];
rz(1.5003381) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7805444) q[0];
sx q[0];
rz(-2.2760976) q[0];
sx q[0];
rz(-1.7678827) q[0];
rz(-0.79958493) q[1];
sx q[1];
rz(-1.0069964) q[1];
sx q[1];
rz(1.132157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06057418) q[0];
sx q[0];
rz(-1.7857528) q[0];
sx q[0];
rz(0.092553986) q[0];
rz(-pi) q[1];
rz(-2.8707809) q[2];
sx q[2];
rz(-3.0584444) q[2];
sx q[2];
rz(1.5025592) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2136176) q[1];
sx q[1];
rz(-2.2384751) q[1];
sx q[1];
rz(-1.4721406) q[1];
x q[2];
rz(1.6145094) q[3];
sx q[3];
rz(-0.96230405) q[3];
sx q[3];
rz(2.9491346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.75480294) q[2];
sx q[2];
rz(-0.58480442) q[2];
sx q[2];
rz(1.3873788) q[2];
rz(-2.7925708) q[3];
sx q[3];
rz(-1.6882221) q[3];
sx q[3];
rz(0.52946985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3997407) q[0];
sx q[0];
rz(-2.1735503) q[0];
sx q[0];
rz(0.35811785) q[0];
rz(1.0707567) q[1];
sx q[1];
rz(-0.64859575) q[1];
sx q[1];
rz(1.4395641) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4442397) q[0];
sx q[0];
rz(-1.8336189) q[0];
sx q[0];
rz(0.79189827) q[0];
rz(-0.87936833) q[2];
sx q[2];
rz(-1.9676529) q[2];
sx q[2];
rz(-0.59890998) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6144952) q[1];
sx q[1];
rz(-2.1367837) q[1];
sx q[1];
rz(2.4742592) q[1];
x q[2];
rz(1.8240806) q[3];
sx q[3];
rz(-0.82013762) q[3];
sx q[3];
rz(1.9092321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4293356) q[2];
sx q[2];
rz(-1.9108994) q[2];
sx q[2];
rz(2.5168391) q[2];
rz(-1.6338927) q[3];
sx q[3];
rz(-2.3816536) q[3];
sx q[3];
rz(1.1225351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.7655012) q[0];
sx q[0];
rz(-2.2608345) q[0];
sx q[0];
rz(1.9951903) q[0];
rz(1.7064077) q[1];
sx q[1];
rz(-2.621666) q[1];
sx q[1];
rz(-0.82040876) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7026742) q[0];
sx q[0];
rz(-1.6587388) q[0];
sx q[0];
rz(-1.8431208) q[0];
x q[1];
rz(1.3455639) q[2];
sx q[2];
rz(-0.67906717) q[2];
sx q[2];
rz(-1.0240384) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31215224) q[1];
sx q[1];
rz(-1.7746801) q[1];
sx q[1];
rz(1.3711434) q[1];
x q[2];
rz(2.5518604) q[3];
sx q[3];
rz(-1.5431649) q[3];
sx q[3];
rz(3.0363135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4121805) q[2];
sx q[2];
rz(-2.086144) q[2];
sx q[2];
rz(0.5385651) q[2];
rz(1.6719079) q[3];
sx q[3];
rz(-1.6939751) q[3];
sx q[3];
rz(-0.058549747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.3013714) q[0];
sx q[0];
rz(-3.0026307) q[0];
sx q[0];
rz(3.050991) q[0];
rz(0.36965707) q[1];
sx q[1];
rz(-1.548111) q[1];
sx q[1];
rz(0.37818092) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9210631) q[0];
sx q[0];
rz(-0.53595966) q[0];
sx q[0];
rz(0.45951636) q[0];
rz(-pi) q[1];
rz(-0.090094968) q[2];
sx q[2];
rz(-1.6125154) q[2];
sx q[2];
rz(-0.99064529) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.52247574) q[1];
sx q[1];
rz(-1.4475087) q[1];
sx q[1];
rz(-0.49501623) q[1];
rz(-pi) q[2];
rz(0.079433283) q[3];
sx q[3];
rz(-1.7609247) q[3];
sx q[3];
rz(-0.058323764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6211264) q[2];
sx q[2];
rz(-1.6660606) q[2];
sx q[2];
rz(-0.10144083) q[2];
rz(1.8488047) q[3];
sx q[3];
rz(-2.3088876) q[3];
sx q[3];
rz(1.4147991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3066278) q[0];
sx q[0];
rz(-1.2766301) q[0];
sx q[0];
rz(-2.2391338) q[0];
rz(-2.9023671) q[1];
sx q[1];
rz(-1.7996412) q[1];
sx q[1];
rz(-0.36144027) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.343924) q[0];
sx q[0];
rz(-1.6822681) q[0];
sx q[0];
rz(2.0693857) q[0];
rz(-pi) q[1];
rz(-0.41255422) q[2];
sx q[2];
rz(-3.0027886) q[2];
sx q[2];
rz(-2.1766162) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4415008) q[1];
sx q[1];
rz(-2.6569286) q[1];
sx q[1];
rz(-1.086471) q[1];
rz(2.7167201) q[3];
sx q[3];
rz(-2.7719404) q[3];
sx q[3];
rz(1.2583789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.07301894) q[2];
sx q[2];
rz(-1.34015) q[2];
sx q[2];
rz(0.87693357) q[2];
rz(-2.1584611) q[3];
sx q[3];
rz(-1.5777595) q[3];
sx q[3];
rz(-2.1985998) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8125732) q[0];
sx q[0];
rz(-0.1952157) q[0];
sx q[0];
rz(-0.83874291) q[0];
rz(2.4328531) q[1];
sx q[1];
rz(-0.63116169) q[1];
sx q[1];
rz(0.90604025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4137461) q[0];
sx q[0];
rz(-2.4165396) q[0];
sx q[0];
rz(3.1097799) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9211414) q[2];
sx q[2];
rz(-1.7736777) q[2];
sx q[2];
rz(1.6914934) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3248277) q[1];
sx q[1];
rz(-0.93597368) q[1];
sx q[1];
rz(-0.88714182) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1791547) q[3];
sx q[3];
rz(-1.3624853) q[3];
sx q[3];
rz(2.359451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9771542) q[2];
sx q[2];
rz(-0.70455569) q[2];
sx q[2];
rz(-1.1158811) q[2];
rz(0.62853938) q[3];
sx q[3];
rz(-2.0597337) q[3];
sx q[3];
rz(-2.5270497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3299265) q[0];
sx q[0];
rz(-2.3217432) q[0];
sx q[0];
rz(2.9803168) q[0];
rz(0.8423841) q[1];
sx q[1];
rz(-1.7120275) q[1];
sx q[1];
rz(2.9715723) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9901491) q[0];
sx q[0];
rz(-3.121576) q[0];
sx q[0];
rz(-0.41949864) q[0];
rz(-1.4030946) q[2];
sx q[2];
rz(-1.1360886) q[2];
sx q[2];
rz(-2.0185883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.051703) q[1];
sx q[1];
rz(-2.0031628) q[1];
sx q[1];
rz(-2.927455) q[1];
x q[2];
rz(-0.028548553) q[3];
sx q[3];
rz(-1.1880837) q[3];
sx q[3];
rz(0.0023502758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1539479) q[2];
sx q[2];
rz(-1.597007) q[2];
sx q[2];
rz(-1.4863996) q[2];
rz(-3.0814643) q[3];
sx q[3];
rz(-2.3190053) q[3];
sx q[3];
rz(-1.7000343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.5104367) q[0];
sx q[0];
rz(-0.031351723) q[0];
sx q[0];
rz(-1.7106868) q[0];
rz(2.778964) q[1];
sx q[1];
rz(-1.5321833) q[1];
sx q[1];
rz(-1.8261725) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5820438) q[0];
sx q[0];
rz(-1.2799731) q[0];
sx q[0];
rz(-1.2469588) q[0];
x q[1];
rz(1.1926629) q[2];
sx q[2];
rz(-1.4571428) q[2];
sx q[2];
rz(3.0494351) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.77389312) q[1];
sx q[1];
rz(-2.4871792) q[1];
sx q[1];
rz(-2.6952637) q[1];
rz(-0.86851991) q[3];
sx q[3];
rz(-1.9316626) q[3];
sx q[3];
rz(-1.3165733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0843087) q[2];
sx q[2];
rz(-1.5949275) q[2];
sx q[2];
rz(0.16051897) q[2];
rz(-2.6594243) q[3];
sx q[3];
rz(-0.67626685) q[3];
sx q[3];
rz(-0.67960656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-3.1228444) q[0];
sx q[0];
rz(-1.7457122) q[0];
sx q[0];
rz(2.2298298) q[0];
rz(-1.1625166) q[1];
sx q[1];
rz(-1.5717506) q[1];
sx q[1];
rz(2.7350978) q[1];
rz(0.087333655) q[2];
sx q[2];
rz(-1.6651911) q[2];
sx q[2];
rz(-1.390425) q[2];
rz(-1.0033458) q[3];
sx q[3];
rz(-0.88450817) q[3];
sx q[3];
rz(-1.3664288) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
