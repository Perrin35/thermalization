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
rz(-3.128669) q[0];
rz(0.68459964) q[1];
sx q[1];
rz(-2.3426988) q[1];
sx q[1];
rz(1.0577143) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81257129) q[0];
sx q[0];
rz(-0.91513915) q[0];
sx q[0];
rz(1.283487) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3637412) q[2];
sx q[2];
rz(-1.9490644) q[2];
sx q[2];
rz(3.0562378) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8243421) q[1];
sx q[1];
rz(-1.2095272) q[1];
sx q[1];
rz(-2.3995705) q[1];
rz(1.9625447) q[3];
sx q[3];
rz(-1.6101675) q[3];
sx q[3];
rz(1.5567428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5554123) q[2];
sx q[2];
rz(-1.5427898) q[2];
sx q[2];
rz(-3.0482698) q[2];
rz(-0.020545067) q[3];
sx q[3];
rz(-13/(16*pi)) q[3];
sx q[3];
rz(1.7790022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-3.0697295) q[0];
sx q[0];
rz(-1.7427895) q[0];
sx q[0];
rz(-2.3139957) q[0];
rz(0.96356511) q[1];
sx q[1];
rz(-1.5255442) q[1];
sx q[1];
rz(0.73659426) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927836) q[0];
sx q[0];
rz(-0.28454706) q[0];
sx q[0];
rz(-1.8440767) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57666333) q[2];
sx q[2];
rz(-2.0915589) q[2];
sx q[2];
rz(-1.3530089) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5399038) q[1];
sx q[1];
rz(-1.772953) q[1];
sx q[1];
rz(-1.8784932) q[1];
rz(-1.5483891) q[3];
sx q[3];
rz(-1.7051892) q[3];
sx q[3];
rz(2.973864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.57614148) q[2];
sx q[2];
rz(-0.57500035) q[2];
sx q[2];
rz(0.74419332) q[2];
rz(-0.60892504) q[3];
sx q[3];
rz(-0.78151339) q[3];
sx q[3];
rz(-1.6412546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3610483) q[0];
sx q[0];
rz(-0.86549509) q[0];
sx q[0];
rz(-1.3737099) q[0];
rz(-0.79958493) q[1];
sx q[1];
rz(-1.0069964) q[1];
sx q[1];
rz(-2.0094357) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47100859) q[0];
sx q[0];
rz(-0.23375227) q[0];
sx q[0];
rz(1.170355) q[0];
rz(3.0614616) q[2];
sx q[2];
rz(-1.5930158) q[2];
sx q[2];
rz(2.9399088) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.76946041) q[1];
sx q[1];
rz(-2.4677708) q[1];
sx q[1];
rz(-3.0173183) q[1];
rz(-pi) q[2];
rz(0.60894062) q[3];
sx q[3];
rz(-1.534933) q[3];
sx q[3];
rz(-1.3533398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3867897) q[2];
sx q[2];
rz(-2.5567882) q[2];
sx q[2];
rz(-1.7542138) q[2];
rz(-0.34902188) q[3];
sx q[3];
rz(-1.6882221) q[3];
sx q[3];
rz(-0.52946985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3997407) q[0];
sx q[0];
rz(-2.1735503) q[0];
sx q[0];
rz(0.35811785) q[0];
rz(2.070836) q[1];
sx q[1];
rz(-2.4929969) q[1];
sx q[1];
rz(-1.7020285) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69735295) q[0];
sx q[0];
rz(-1.3079738) q[0];
sx q[0];
rz(2.3496944) q[0];
rz(-pi) q[1];
rz(-2.6433377) q[2];
sx q[2];
rz(-2.1994805) q[2];
sx q[2];
rz(-1.2815338) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.52709748) q[1];
sx q[1];
rz(-1.004809) q[1];
sx q[1];
rz(-0.66733349) q[1];
rz(1.8240806) q[3];
sx q[3];
rz(-2.321455) q[3];
sx q[3];
rz(-1.9092321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4293356) q[2];
sx q[2];
rz(-1.2306932) q[2];
sx q[2];
rz(-0.62475359) q[2];
rz(-1.5077) q[3];
sx q[3];
rz(-2.3816536) q[3];
sx q[3];
rz(-1.1225351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.7655012) q[0];
sx q[0];
rz(-2.2608345) q[0];
sx q[0];
rz(-1.9951903) q[0];
rz(1.435185) q[1];
sx q[1];
rz(-2.621666) q[1];
sx q[1];
rz(0.82040876) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2489481) q[0];
sx q[0];
rz(-1.8420418) q[0];
sx q[0];
rz(0.091288996) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3455639) q[2];
sx q[2];
rz(-2.4625255) q[2];
sx q[2];
rz(-2.1175543) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0972406) q[1];
sx q[1];
rz(-2.8572081) q[1];
sx q[1];
rz(0.76459717) q[1];
rz(-0.049656258) q[3];
sx q[3];
rz(-0.59030246) q[3];
sx q[3];
rz(1.4242537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.72941214) q[2];
sx q[2];
rz(-1.0554487) q[2];
sx q[2];
rz(2.6030276) q[2];
rz(-1.4696848) q[3];
sx q[3];
rz(-1.4476176) q[3];
sx q[3];
rz(0.058549747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3013714) q[0];
sx q[0];
rz(-3.0026307) q[0];
sx q[0];
rz(0.090601623) q[0];
rz(0.36965707) q[1];
sx q[1];
rz(-1.5934817) q[1];
sx q[1];
rz(-0.37818092) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3017186) q[0];
sx q[0];
rz(-1.0953971) q[0];
sx q[0];
rz(-1.8283707) q[0];
rz(0.090094968) q[2];
sx q[2];
rz(-1.5290773) q[2];
sx q[2];
rz(2.1509474) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6191169) q[1];
sx q[1];
rz(-1.694084) q[1];
sx q[1];
rz(-0.49501623) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0621594) q[3];
sx q[3];
rz(-1.380668) q[3];
sx q[3];
rz(-3.0832689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6211264) q[2];
sx q[2];
rz(-1.4755321) q[2];
sx q[2];
rz(-3.0401518) q[2];
rz(-1.8488047) q[3];
sx q[3];
rz(-0.83270508) q[3];
sx q[3];
rz(-1.7267936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8349649) q[0];
sx q[0];
rz(-1.8649626) q[0];
sx q[0];
rz(-2.2391338) q[0];
rz(0.2392256) q[1];
sx q[1];
rz(-1.7996412) q[1];
sx q[1];
rz(2.7801524) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.343924) q[0];
sx q[0];
rz(-1.4593246) q[0];
sx q[0];
rz(1.072207) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5148411) q[2];
sx q[2];
rz(-1.6978886) q[2];
sx q[2];
rz(-2.5927319) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8351961) q[1];
sx q[1];
rz(-1.3521242) q[1];
sx q[1];
rz(-1.1347215) q[1];
rz(-pi) q[2];
rz(-2.802235) q[3];
sx q[3];
rz(-1.7202783) q[3];
sx q[3];
rz(0.71163346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0685737) q[2];
sx q[2];
rz(-1.34015) q[2];
sx q[2];
rz(-0.87693357) q[2];
rz(2.1584611) q[3];
sx q[3];
rz(-1.5777595) q[3];
sx q[3];
rz(2.1985998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8125732) q[0];
sx q[0];
rz(-2.946377) q[0];
sx q[0];
rz(0.83874291) q[0];
rz(2.4328531) q[1];
sx q[1];
rz(-0.63116169) q[1];
sx q[1];
rz(0.90604025) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7278465) q[0];
sx q[0];
rz(-2.4165396) q[0];
sx q[0];
rz(0.031812761) q[0];
rz(-1.2204513) q[2];
sx q[2];
rz(-1.7736777) q[2];
sx q[2];
rz(-1.4500993) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4455394) q[1];
sx q[1];
rz(-2.1043244) q[1];
sx q[1];
rz(2.3818156) q[1];
x q[2];
rz(2.8894807) q[3];
sx q[3];
rz(-0.97741717) q[3];
sx q[3];
rz(-0.93170792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.16443843) q[2];
sx q[2];
rz(-2.437037) q[2];
sx q[2];
rz(2.0257115) q[2];
rz(-0.62853938) q[3];
sx q[3];
rz(-2.0597337) q[3];
sx q[3];
rz(-0.61454296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3299265) q[0];
sx q[0];
rz(-0.8198494) q[0];
sx q[0];
rz(0.16127583) q[0];
rz(2.2992086) q[1];
sx q[1];
rz(-1.4295652) q[1];
sx q[1];
rz(-0.17002034) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5710167) q[0];
sx q[0];
rz(-1.5890772) q[0];
sx q[0];
rz(-1.5789501) q[0];
rz(0.34509322) q[2];
sx q[2];
rz(-2.6775914) q[2];
sx q[2];
rz(0.74079266) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0898897) q[1];
sx q[1];
rz(-2.0031628) q[1];
sx q[1];
rz(0.21413762) q[1];
rz(-pi) q[2];
rz(-1.6415855) q[3];
sx q[3];
rz(-0.38372358) q[3];
sx q[3];
rz(3.0629223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9876447) q[2];
sx q[2];
rz(-1.597007) q[2];
sx q[2];
rz(-1.4863996) q[2];
rz(-0.060128309) q[3];
sx q[3];
rz(-2.3190053) q[3];
sx q[3];
rz(-1.4415584) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63115591) q[0];
sx q[0];
rz(-3.1102409) q[0];
sx q[0];
rz(-1.4309058) q[0];
rz(0.36262861) q[1];
sx q[1];
rz(-1.6094094) q[1];
sx q[1];
rz(1.3154202) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5595488) q[0];
sx q[0];
rz(-1.8616195) q[0];
sx q[0];
rz(-1.2469588) q[0];
x q[1];
rz(-1.1926629) q[2];
sx q[2];
rz(-1.6844498) q[2];
sx q[2];
rz(-0.092157539) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3676995) q[1];
sx q[1];
rz(-2.4871792) q[1];
sx q[1];
rz(-2.6952637) q[1];
x q[2];
rz(-2.2730727) q[3];
sx q[3];
rz(-1.2099301) q[3];
sx q[3];
rz(-1.3165733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.057283904) q[2];
sx q[2];
rz(-1.5466651) q[2];
sx q[2];
rz(-2.9810737) q[2];
rz(-2.6594243) q[3];
sx q[3];
rz(-0.67626685) q[3];
sx q[3];
rz(2.4619861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1228444) q[0];
sx q[0];
rz(-1.7457122) q[0];
sx q[0];
rz(2.2298298) q[0];
rz(1.979076) q[1];
sx q[1];
rz(-1.5717506) q[1];
sx q[1];
rz(2.7350978) q[1];
rz(0.82635098) q[2];
sx q[2];
rz(-3.0130825) q[2];
sx q[2];
rz(2.4997406) q[2];
rz(-2.1382469) q[3];
sx q[3];
rz(-2.2570845) q[3];
sx q[3];
rz(1.7751638) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
