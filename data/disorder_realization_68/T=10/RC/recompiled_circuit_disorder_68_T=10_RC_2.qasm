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
rz(-0.36748537) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4700692) q[0];
sx q[0];
rz(-2.1052261) q[0];
sx q[0];
rz(0.1202017) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3221402) q[2];
sx q[2];
rz(-2.6373632) q[2];
sx q[2];
rz(2.3766975) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6403113) q[1];
sx q[1];
rz(-2.8456563) q[1];
sx q[1];
rz(-2.179115) q[1];
rz(2.8005881) q[3];
sx q[3];
rz(-1.4031646) q[3];
sx q[3];
rz(-2.117702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.964103) q[2];
sx q[2];
rz(-2.6314645) q[2];
sx q[2];
rz(-2.5906079) q[2];
rz(1.3059113) q[3];
sx q[3];
rz(-1.4923613) q[3];
sx q[3];
rz(1.3163542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47857639) q[0];
sx q[0];
rz(-1.0000279) q[0];
sx q[0];
rz(2.6696894) q[0];
rz(2.7117803) q[1];
sx q[1];
rz(-1.8919573) q[1];
sx q[1];
rz(2.205251) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72915709) q[0];
sx q[0];
rz(-1.658174) q[0];
sx q[0];
rz(-0.23075128) q[0];
rz(-0.78511946) q[2];
sx q[2];
rz(-1.2650507) q[2];
sx q[2];
rz(2.608125) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1658926) q[1];
sx q[1];
rz(-0.70043889) q[1];
sx q[1];
rz(-0.16209929) q[1];
rz(-1.068088) q[3];
sx q[3];
rz(-1.5912676) q[3];
sx q[3];
rz(0.81077829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.77461809) q[2];
sx q[2];
rz(-2.8145511) q[2];
sx q[2];
rz(2.7152087) q[2];
rz(-1.9042227) q[3];
sx q[3];
rz(-0.62785134) q[3];
sx q[3];
rz(3.1085076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.8957829) q[0];
sx q[0];
rz(-1.3170467) q[0];
sx q[0];
rz(0.93908969) q[0];
rz(2.242873) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(-0.59392196) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9575189) q[0];
sx q[0];
rz(-2.6960417) q[0];
sx q[0];
rz(0.81934388) q[0];
x q[1];
rz(0.8823231) q[2];
sx q[2];
rz(-1.6154628) q[2];
sx q[2];
rz(-0.67827144) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3214026) q[1];
sx q[1];
rz(-0.47649511) q[1];
sx q[1];
rz(2.0501775) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61641683) q[3];
sx q[3];
rz(-0.55509242) q[3];
sx q[3];
rz(-2.4544231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5014191) q[2];
sx q[2];
rz(-2.4121425) q[2];
sx q[2];
rz(-1.7017986) q[2];
rz(-2.7539608) q[3];
sx q[3];
rz(-1.6250316) q[3];
sx q[3];
rz(-0.38813996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7664292) q[0];
sx q[0];
rz(-1.5513993) q[0];
sx q[0];
rz(0.50278062) q[0];
rz(-2.373383) q[1];
sx q[1];
rz(-2.6380824) q[1];
sx q[1];
rz(2.3847413) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4810281) q[0];
sx q[0];
rz(-1.2111944) q[0];
sx q[0];
rz(1.9268376) q[0];
x q[1];
rz(1.1050622) q[2];
sx q[2];
rz(-1.5591991) q[2];
sx q[2];
rz(2.2955745) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6824324) q[1];
sx q[1];
rz(-1.2130514) q[1];
sx q[1];
rz(-2.1898502) q[1];
x q[2];
rz(0.54550708) q[3];
sx q[3];
rz(-0.81199284) q[3];
sx q[3];
rz(0.10520392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.42671529) q[2];
sx q[2];
rz(-1.2299512) q[2];
sx q[2];
rz(1.4871917) q[2];
rz(-0.58250827) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(0.55707651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8476167) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(-2.3838682) q[0];
rz(1.853653) q[1];
sx q[1];
rz(-2.2133591) q[1];
sx q[1];
rz(2.0910738) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3462853) q[0];
sx q[0];
rz(-1.8827569) q[0];
sx q[0];
rz(0.26766582) q[0];
x q[1];
rz(0.96111416) q[2];
sx q[2];
rz(-2.4498307) q[2];
sx q[2];
rz(1.8209396) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.95549059) q[1];
sx q[1];
rz(-2.0932066) q[1];
sx q[1];
rz(0.91764692) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53253048) q[3];
sx q[3];
rz(-2.2610287) q[3];
sx q[3];
rz(-1.1158451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.22333764) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(0.63344947) q[2];
rz(-1.194362) q[3];
sx q[3];
rz(-1.6703689) q[3];
sx q[3];
rz(0.71715322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.441992) q[0];
sx q[0];
rz(-0.49015912) q[0];
sx q[0];
rz(-2.8884086) q[0];
rz(1.5340012) q[1];
sx q[1];
rz(-1.7065159) q[1];
sx q[1];
rz(-1.6794499) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7397241) q[0];
sx q[0];
rz(-1.7812294) q[0];
sx q[0];
rz(1.6519288) q[0];
x q[1];
rz(2.9003733) q[2];
sx q[2];
rz(-0.93572817) q[2];
sx q[2];
rz(-1.5347753) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6219382) q[1];
sx q[1];
rz(-0.42236537) q[1];
sx q[1];
rz(-2.7154891) q[1];
rz(-pi) q[2];
rz(-2.3092689) q[3];
sx q[3];
rz(-0.451085) q[3];
sx q[3];
rz(3.0576599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9399461) q[2];
sx q[2];
rz(-0.74735171) q[2];
sx q[2];
rz(0.80491006) q[2];
rz(1.1770052) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(-2.0578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(-0.2729623) q[0];
sx q[0];
rz(-2.0721764) q[0];
sx q[0];
rz(-0.92765635) q[0];
rz(1.0246798) q[1];
sx q[1];
rz(-1.506348) q[1];
sx q[1];
rz(1.0120846) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18400684) q[0];
sx q[0];
rz(-1.1054966) q[0];
sx q[0];
rz(2.1561949) q[0];
rz(-1.4321248) q[2];
sx q[2];
rz(-1.1853293) q[2];
sx q[2];
rz(-0.0027545714) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.63230292) q[1];
sx q[1];
rz(-1.0273233) q[1];
sx q[1];
rz(1.6093045) q[1];
rz(-0.13109644) q[3];
sx q[3];
rz(-2.5790865) q[3];
sx q[3];
rz(-1.8545811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6107789) q[2];
sx q[2];
rz(-1.4799708) q[2];
sx q[2];
rz(-2.7116595) q[2];
rz(1.0144462) q[3];
sx q[3];
rz(-0.40922624) q[3];
sx q[3];
rz(2.608192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5291418) q[0];
sx q[0];
rz(-0.93944678) q[0];
sx q[0];
rz(-2.9274143) q[0];
rz(2.0902436) q[1];
sx q[1];
rz(-0.21251692) q[1];
sx q[1];
rz(-0.28373757) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8026233) q[0];
sx q[0];
rz(-1.5366652) q[0];
sx q[0];
rz(0.3011093) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19212888) q[2];
sx q[2];
rz(-1.8953952) q[2];
sx q[2];
rz(-1.0629551) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6730496) q[1];
sx q[1];
rz(-1.0769516) q[1];
sx q[1];
rz(1.3480575) q[1];
rz(-pi) q[2];
rz(-0.92054263) q[3];
sx q[3];
rz(-2.1170108) q[3];
sx q[3];
rz(-1.9073245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6909137) q[2];
sx q[2];
rz(-0.51270715) q[2];
sx q[2];
rz(1.696375) q[2];
rz(1.5971659) q[3];
sx q[3];
rz(-1.701136) q[3];
sx q[3];
rz(0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.504869) q[0];
sx q[0];
rz(-3.0915785) q[0];
sx q[0];
rz(-0.069256393) q[0];
rz(1.6537369) q[1];
sx q[1];
rz(-1.8585049) q[1];
sx q[1];
rz(-1.5690631) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3587787) q[0];
sx q[0];
rz(-0.94952119) q[0];
sx q[0];
rz(0.57399477) q[0];
rz(2.9022129) q[2];
sx q[2];
rz(-0.33472543) q[2];
sx q[2];
rz(0.3604381) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3410586) q[1];
sx q[1];
rz(-1.8039928) q[1];
sx q[1];
rz(-1.6569767) q[1];
rz(-0.2089573) q[3];
sx q[3];
rz(-2.2088802) q[3];
sx q[3];
rz(1.0554505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9562324) q[2];
sx q[2];
rz(-2.9113443) q[2];
sx q[2];
rz(0.13988477) q[2];
rz(-2.774003) q[3];
sx q[3];
rz(-1.9544173) q[3];
sx q[3];
rz(2.1504413) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763828) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(0.64176732) q[0];
rz(1.9104674) q[1];
sx q[1];
rz(-1.9893913) q[1];
sx q[1];
rz(-0.26783255) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56775996) q[0];
sx q[0];
rz(-1.1288252) q[0];
sx q[0];
rz(-2.3964336) q[0];
x q[1];
rz(-0.97271131) q[2];
sx q[2];
rz(-0.22062606) q[2];
sx q[2];
rz(1.0704744) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8903058) q[1];
sx q[1];
rz(-2.7000513) q[1];
sx q[1];
rz(2.4114386) q[1];
rz(-pi) q[2];
rz(2.855741) q[3];
sx q[3];
rz(-2.0945858) q[3];
sx q[3];
rz(1.3956192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.81007593) q[2];
sx q[2];
rz(-1.6167567) q[2];
sx q[2];
rz(-1.6646741) q[2];
rz(-0.26633513) q[3];
sx q[3];
rz(-2.895152) q[3];
sx q[3];
rz(0.5464856) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01263604) q[0];
sx q[0];
rz(-0.92100443) q[0];
sx q[0];
rz(-2.2367649) q[0];
rz(2.3616882) q[1];
sx q[1];
rz(-0.48702469) q[1];
sx q[1];
rz(-1.3866966) q[1];
rz(-0.98942479) q[2];
sx q[2];
rz(-2.1944254) q[2];
sx q[2];
rz(2.2133322) q[2];
rz(-0.6775425) q[3];
sx q[3];
rz(-1.7384221) q[3];
sx q[3];
rz(1.3330028) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
