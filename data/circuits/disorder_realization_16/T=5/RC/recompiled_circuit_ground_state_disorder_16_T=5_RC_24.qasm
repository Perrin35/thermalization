OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0082173) q[0];
sx q[0];
rz(-1.2674588) q[0];
sx q[0];
rz(-0.01292364) q[0];
rz(0.68459964) q[1];
sx q[1];
rz(-2.3426988) q[1];
sx q[1];
rz(1.0577143) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36127451) q[0];
sx q[0];
rz(-2.4343505) q[0];
sx q[0];
rz(-2.7886084) q[0];
rz(-pi) q[1];
rz(2.3637412) q[2];
sx q[2];
rz(-1.9490644) q[2];
sx q[2];
rz(0.08535484) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.94016852) q[1];
sx q[1];
rz(-2.2551309) q[1];
sx q[1];
rz(-1.0971054) q[1];
rz(-pi) q[2];
rz(-1.9625447) q[3];
sx q[3];
rz(-1.5314252) q[3];
sx q[3];
rz(-1.5848499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.58618033) q[2];
sx q[2];
rz(-1.5427898) q[2];
sx q[2];
rz(3.0482698) q[2];
rz(3.1210476) q[3];
sx q[3];
rz(-2.8829657) q[3];
sx q[3];
rz(-1.7790022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0697295) q[0];
sx q[0];
rz(-1.3988031) q[0];
sx q[0];
rz(-2.3139957) q[0];
rz(-0.96356511) q[1];
sx q[1];
rz(-1.5255442) q[1];
sx q[1];
rz(2.4049984) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4086491) q[0];
sx q[0];
rz(-1.8445065) q[0];
sx q[0];
rz(-3.0628171) q[0];
rz(-pi) q[1];
rz(2.170855) q[2];
sx q[2];
rz(-2.063437) q[2];
sx q[2];
rz(-0.095183177) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40562661) q[1];
sx q[1];
rz(-0.36639226) q[1];
sx q[1];
rz(-0.97586164) q[1];
rz(-pi) q[2];
rz(-3.0071665) q[3];
sx q[3];
rz(-1.5485912) q[3];
sx q[3];
rz(-1.7355222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5654512) q[2];
sx q[2];
rz(-0.57500035) q[2];
sx q[2];
rz(-2.3973993) q[2];
rz(-0.60892504) q[3];
sx q[3];
rz(-0.78151339) q[3];
sx q[3];
rz(-1.6412546) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3610483) q[0];
sx q[0];
rz(-0.86549509) q[0];
sx q[0];
rz(-1.7678827) q[0];
rz(-2.3420077) q[1];
sx q[1];
rz(-2.1345963) q[1];
sx q[1];
rz(-2.0094357) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06057418) q[0];
sx q[0];
rz(-1.3558398) q[0];
sx q[0];
rz(-0.092553986) q[0];
rz(0.08013101) q[2];
sx q[2];
rz(-1.5930158) q[2];
sx q[2];
rz(0.20168389) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3721322) q[1];
sx q[1];
rz(-0.67382183) q[1];
sx q[1];
rz(0.12427434) q[1];
rz(-pi) q[2];
x q[2];
rz(2.532652) q[3];
sx q[3];
rz(-1.6066597) q[3];
sx q[3];
rz(-1.3533398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.75480294) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.741852) q[0];
sx q[0];
rz(-2.1735503) q[0];
sx q[0];
rz(-0.35811785) q[0];
rz(-2.070836) q[1];
sx q[1];
rz(-0.64859575) q[1];
sx q[1];
rz(-1.7020285) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5191906) q[0];
sx q[0];
rz(-0.82525142) q[0];
sx q[0];
rz(-2.78016) q[0];
rz(-pi) q[1];
rz(2.1522572) q[2];
sx q[2];
rz(-0.78063595) q[2];
sx q[2];
rz(-0.5350185) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4434699) q[1];
sx q[1];
rz(-1.021153) q[1];
sx q[1];
rz(-0.89068954) q[1];
rz(-pi) q[2];
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
rz(-1.7122571) q[2];
sx q[2];
rz(-1.2306932) q[2];
sx q[2];
rz(2.5168391) q[2];
rz(-1.5077) q[3];
sx q[3];
rz(-0.75993901) q[3];
sx q[3];
rz(1.1225351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7655012) q[0];
sx q[0];
rz(-2.2608345) q[0];
sx q[0];
rz(1.1464024) q[0];
rz(1.435185) q[1];
sx q[1];
rz(-0.51992661) q[1];
sx q[1];
rz(-0.82040876) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89264458) q[0];
sx q[0];
rz(-1.2995509) q[0];
sx q[0];
rz(-3.0503037) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7960288) q[2];
sx q[2];
rz(-0.67906717) q[2];
sx q[2];
rz(1.0240384) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0972406) q[1];
sx q[1];
rz(-2.8572081) q[1];
sx q[1];
rz(2.3769955) q[1];
x q[2];
rz(1.6040398) q[3];
sx q[3];
rz(-2.1602732) q[3];
sx q[3];
rz(-1.4840028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72941214) q[2];
sx q[2];
rz(-2.086144) q[2];
sx q[2];
rz(0.5385651) q[2];
rz(-1.6719079) q[3];
sx q[3];
rz(-1.4476176) q[3];
sx q[3];
rz(3.0830429) q[3];
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
rz(2.3013714) q[0];
sx q[0];
rz(-0.138962) q[0];
sx q[0];
rz(0.090601623) q[0];
rz(2.7719356) q[1];
sx q[1];
rz(-1.548111) q[1];
sx q[1];
rz(-0.37818092) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3017186) q[0];
sx q[0];
rz(-2.0461956) q[0];
sx q[0];
rz(1.8283707) q[0];
rz(-0.090094968) q[2];
sx q[2];
rz(-1.6125154) q[2];
sx q[2];
rz(-0.99064529) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52247574) q[1];
sx q[1];
rz(-1.4475087) q[1];
sx q[1];
rz(-2.6465764) q[1];
rz(-pi) q[2];
rz(3.0621594) q[3];
sx q[3];
rz(-1.380668) q[3];
sx q[3];
rz(-0.058323764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6211264) q[2];
sx q[2];
rz(-1.4755321) q[2];
sx q[2];
rz(-3.0401518) q[2];
rz(-1.2927879) q[3];
sx q[3];
rz(-0.83270508) q[3];
sx q[3];
rz(-1.4147991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3066278) q[0];
sx q[0];
rz(-1.8649626) q[0];
sx q[0];
rz(2.2391338) q[0];
rz(0.2392256) q[1];
sx q[1];
rz(-1.7996412) q[1];
sx q[1];
rz(-0.36144027) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1162565) q[0];
sx q[0];
rz(-0.5098719) q[0];
sx q[0];
rz(-1.340853) q[0];
rz(-0.41255422) q[2];
sx q[2];
rz(-3.0027886) q[2];
sx q[2];
rz(-2.1766162) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.16365227) q[1];
sx q[1];
rz(-1.1457902) q[1];
sx q[1];
rz(0.24042701) q[1];
rz(-pi) q[2];
rz(-2.802235) q[3];
sx q[3];
rz(-1.7202783) q[3];
sx q[3];
rz(0.71163346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.07301894) q[2];
sx q[2];
rz(-1.8014427) q[2];
sx q[2];
rz(0.87693357) q[2];
rz(2.1584611) q[3];
sx q[3];
rz(-1.5777595) q[3];
sx q[3];
rz(-0.94299281) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32901949) q[0];
sx q[0];
rz(-0.1952157) q[0];
sx q[0];
rz(0.83874291) q[0];
rz(-0.70873952) q[1];
sx q[1];
rz(-0.63116169) q[1];
sx q[1];
rz(-2.2355524) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6853537) q[0];
sx q[0];
rz(-2.2954012) q[0];
sx q[0];
rz(1.5426209) q[0];
rz(-1.0308517) q[2];
sx q[2];
rz(-2.738852) q[2];
sx q[2];
rz(2.7583721) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4455394) q[1];
sx q[1];
rz(-2.1043244) q[1];
sx q[1];
rz(2.3818156) q[1];
x q[2];
rz(0.25211199) q[3];
sx q[3];
rz(-0.97741717) q[3];
sx q[3];
rz(-2.2098847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9771542) q[2];
sx q[2];
rz(-2.437037) q[2];
sx q[2];
rz(-2.0257115) q[2];
rz(0.62853938) q[3];
sx q[3];
rz(-1.081859) q[3];
sx q[3];
rz(2.5270497) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3299265) q[0];
sx q[0];
rz(-2.3217432) q[0];
sx q[0];
rz(0.16127583) q[0];
rz(-2.2992086) q[1];
sx q[1];
rz(-1.7120275) q[1];
sx q[1];
rz(2.9715723) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(7.13037e-05) q[0];
sx q[0];
rz(-1.562644) q[0];
sx q[0];
rz(-0.018281451) q[0];
rz(1.7384981) q[2];
sx q[2];
rz(-2.005504) q[2];
sx q[2];
rz(2.0185883) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5698118) q[1];
sx q[1];
rz(-1.3766292) q[1];
sx q[1];
rz(1.1295736) q[1];
rz(-pi) q[2];
rz(-0.028548553) q[3];
sx q[3];
rz(-1.953509) q[3];
sx q[3];
rz(3.1392424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9876447) q[2];
sx q[2];
rz(-1.5445856) q[2];
sx q[2];
rz(-1.4863996) q[2];
rz(3.0814643) q[3];
sx q[3];
rz(-0.82258737) q[3];
sx q[3];
rz(1.4415584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5104367) q[0];
sx q[0];
rz(-3.1102409) q[0];
sx q[0];
rz(1.7106868) q[0];
rz(-0.36262861) q[1];
sx q[1];
rz(-1.6094094) q[1];
sx q[1];
rz(-1.3154202) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5820438) q[0];
sx q[0];
rz(-1.2799731) q[0];
sx q[0];
rz(-1.2469588) q[0];
rz(-pi) q[1];
rz(-1.8706546) q[2];
sx q[2];
rz(-2.747537) q[2];
sx q[2];
rz(1.9410417) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.23116701) q[1];
sx q[1];
rz(-2.1520414) q[1];
sx q[1];
rz(-1.2509996) q[1];
rz(-pi) q[2];
rz(-0.86851991) q[3];
sx q[3];
rz(-1.2099301) q[3];
sx q[3];
rz(-1.8250193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.057283904) q[2];
sx q[2];
rz(-1.5466651) q[2];
sx q[2];
rz(2.9810737) q[2];
rz(-0.48216835) q[3];
sx q[3];
rz(-2.4653258) q[3];
sx q[3];
rz(2.4619861) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(2.3152417) q[2];
sx q[2];
rz(-0.12851014) q[2];
sx q[2];
rz(-0.64185206) q[2];
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
