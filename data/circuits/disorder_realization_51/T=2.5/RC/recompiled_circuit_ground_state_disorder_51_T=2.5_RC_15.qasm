OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2188323) q[0];
sx q[0];
rz(-0.063475944) q[0];
sx q[0];
rz(1.7712964) q[0];
rz(-2.3308347) q[1];
sx q[1];
rz(-1.6914565) q[1];
sx q[1];
rz(2.9210747) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4654601) q[0];
sx q[0];
rz(-1.5031843) q[0];
sx q[0];
rz(-0.83444447) q[0];
rz(-pi) q[1];
rz(-2.8120997) q[2];
sx q[2];
rz(-1.1348083) q[2];
sx q[2];
rz(1.5832251) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0956362) q[1];
sx q[1];
rz(-1.7929309) q[1];
sx q[1];
rz(-1.5595116) q[1];
x q[2];
rz(-1.0706717) q[3];
sx q[3];
rz(-2.7410738) q[3];
sx q[3];
rz(0.60067117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0567131) q[2];
sx q[2];
rz(-0.006338174) q[2];
sx q[2];
rz(-0.81981266) q[2];
rz(-0.020717185) q[3];
sx q[3];
rz(-0.31532225) q[3];
sx q[3];
rz(2.0899541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3725975) q[0];
sx q[0];
rz(-2.9561434) q[0];
sx q[0];
rz(0.73082596) q[0];
rz(-1.7164879) q[1];
sx q[1];
rz(-2.4162879) q[1];
sx q[1];
rz(-0.46754974) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9502724) q[0];
sx q[0];
rz(-1.5400271) q[0];
sx q[0];
rz(-1.4683718) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5376925) q[2];
sx q[2];
rz(-1.1611801) q[2];
sx q[2];
rz(-2.8996301) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4380355) q[1];
sx q[1];
rz(-2.272637) q[1];
sx q[1];
rz(-1.3753533) q[1];
rz(-pi) q[2];
rz(-2.8465956) q[3];
sx q[3];
rz(-2.2824691) q[3];
sx q[3];
rz(-2.0509348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9759489) q[2];
sx q[2];
rz(-0.011198137) q[2];
sx q[2];
rz(-0.33640081) q[2];
rz(0.30763704) q[3];
sx q[3];
rz(-3.1308789) q[3];
sx q[3];
rz(0.78905869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(0.12863185) q[0];
sx q[0];
rz(-0.19879453) q[0];
sx q[0];
rz(-3.0096753) q[0];
rz(1.4384653) q[1];
sx q[1];
rz(-1.4142282) q[1];
sx q[1];
rz(1.6579113) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8298873) q[0];
sx q[0];
rz(-0.97011203) q[0];
sx q[0];
rz(1.9102824) q[0];
x q[1];
rz(-2.5222579) q[2];
sx q[2];
rz(-0.029994596) q[2];
sx q[2];
rz(-0.66823792) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.10569948) q[1];
sx q[1];
rz(-1.4669682) q[1];
sx q[1];
rz(0.052171252) q[1];
x q[2];
rz(2.4529732) q[3];
sx q[3];
rz(-0.11121777) q[3];
sx q[3];
rz(-2.2764549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.57054532) q[2];
sx q[2];
rz(-1.3344301) q[2];
sx q[2];
rz(1.6837616) q[2];
rz(-0.7782065) q[3];
sx q[3];
rz(-0.0056548803) q[3];
sx q[3];
rz(-2.9290504) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0197765) q[0];
sx q[0];
rz(-1.7541405) q[0];
sx q[0];
rz(-2.8507267) q[0];
rz(-1.5485171) q[1];
sx q[1];
rz(-0.80268186) q[1];
sx q[1];
rz(-0.026550857) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3817859) q[0];
sx q[0];
rz(-1.9142173) q[0];
sx q[0];
rz(1.9476101) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5174116) q[2];
sx q[2];
rz(-2.0211271) q[2];
sx q[2];
rz(0.93102294) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.091152126) q[1];
sx q[1];
rz(-0.9602957) q[1];
sx q[1];
rz(-2.1679085) q[1];
rz(-pi) q[2];
rz(-3.13843) q[3];
sx q[3];
rz(-1.5517276) q[3];
sx q[3];
rz(0.17584383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.87963858) q[2];
sx q[2];
rz(-3.1231572) q[2];
sx q[2];
rz(2.7035942) q[2];
rz(2.0336464) q[3];
sx q[3];
rz(-3.1236533) q[3];
sx q[3];
rz(1.7855135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9172025) q[0];
sx q[0];
rz(-0.49773911) q[0];
sx q[0];
rz(-2.9955731) q[0];
rz(0.13305013) q[1];
sx q[1];
rz(-1.0597884) q[1];
sx q[1];
rz(-0.45566794) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5815703) q[0];
sx q[0];
rz(-3.0455493) q[0];
sx q[0];
rz(0.21012975) q[0];
rz(-2.7002522) q[2];
sx q[2];
rz(-2.5020863) q[2];
sx q[2];
rz(2.9616982) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.86381493) q[1];
sx q[1];
rz(-1.3425273) q[1];
sx q[1];
rz(-0.14240188) q[1];
rz(-pi) q[2];
x q[2];
rz(0.066448575) q[3];
sx q[3];
rz(-0.22804582) q[3];
sx q[3];
rz(1.4888637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2666152) q[2];
sx q[2];
rz(-0.019860331) q[2];
sx q[2];
rz(1.88545) q[2];
rz(2.6202294) q[3];
sx q[3];
rz(-0.25786906) q[3];
sx q[3];
rz(-2.7969587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94619036) q[0];
sx q[0];
rz(-0.38637105) q[0];
sx q[0];
rz(-0.56889164) q[0];
rz(-0.16828123) q[1];
sx q[1];
rz(-1.5852837) q[1];
sx q[1];
rz(0.010244244) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5880347) q[0];
sx q[0];
rz(-1.7419683) q[0];
sx q[0];
rz(1.6200911) q[0];
x q[1];
rz(1.3902789) q[2];
sx q[2];
rz(-1.5897993) q[2];
sx q[2];
rz(-0.190616) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.78624397) q[1];
sx q[1];
rz(-0.054844347) q[1];
sx q[1];
rz(-0.88990239) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5074309) q[3];
sx q[3];
rz(-1.5463357) q[3];
sx q[3];
rz(1.1220049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9804618) q[2];
sx q[2];
rz(-2.9176596) q[2];
sx q[2];
rz(-1.0978318) q[2];
rz(2.791642) q[3];
sx q[3];
rz(-1.8643458) q[3];
sx q[3];
rz(0.85489517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8769281) q[0];
sx q[0];
rz(-0.16328891) q[0];
sx q[0];
rz(0.40633416) q[0];
rz(1.3111275) q[1];
sx q[1];
rz(-2.6268112) q[1];
sx q[1];
rz(-3.0313671) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12719181) q[0];
sx q[0];
rz(-1.4716291) q[0];
sx q[0];
rz(-1.5983461) q[0];
rz(-pi) q[1];
rz(3.1378205) q[2];
sx q[2];
rz(-1.5654148) q[2];
sx q[2];
rz(0.37056915) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.71512039) q[1];
sx q[1];
rz(-2.7139152) q[1];
sx q[1];
rz(-1.5918455) q[1];
rz(-pi) q[2];
rz(-0.87894999) q[3];
sx q[3];
rz(-0.96262299) q[3];
sx q[3];
rz(1.5076306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.280764) q[2];
sx q[2];
rz(-3.132756) q[2];
sx q[2];
rz(-0.78738085) q[2];
rz(-2.8619838) q[3];
sx q[3];
rz(-0.044737261) q[3];
sx q[3];
rz(-0.71981049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1396007) q[0];
sx q[0];
rz(-2.3878492) q[0];
sx q[0];
rz(1.4308223) q[0];
rz(0.17499533) q[1];
sx q[1];
rz(-0.7376968) q[1];
sx q[1];
rz(-1.4281248) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5663877) q[0];
sx q[0];
rz(-0.022179929) q[0];
sx q[0];
rz(-1.7087858) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1375132) q[2];
sx q[2];
rz(-1.578184) q[2];
sx q[2];
rz(-1.0122062) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97615759) q[1];
sx q[1];
rz(-1.7320898) q[1];
sx q[1];
rz(1.2978202) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.346502) q[3];
sx q[3];
rz(-1.845775) q[3];
sx q[3];
rz(-2.0338273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.78508198) q[2];
sx q[2];
rz(-2.3729237) q[2];
sx q[2];
rz(-1.6016426) q[2];
rz(-0.29940638) q[3];
sx q[3];
rz(-1.6573903) q[3];
sx q[3];
rz(-0.33777344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47141075) q[0];
sx q[0];
rz(-3.1306559) q[0];
sx q[0];
rz(-0.48728824) q[0];
rz(1.4393073) q[1];
sx q[1];
rz(-2.3479562) q[1];
sx q[1];
rz(-0.38584858) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7283994) q[0];
sx q[0];
rz(-1.5110713) q[0];
sx q[0];
rz(-0.3932391) q[0];
rz(1.3273238) q[2];
sx q[2];
rz(-0.43691942) q[2];
sx q[2];
rz(-1.5658558) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6035406) q[1];
sx q[1];
rz(-1.5376989) q[1];
sx q[1];
rz(-2.4500928) q[1];
rz(-pi) q[2];
rz(-3.039592) q[3];
sx q[3];
rz(-2.8047562) q[3];
sx q[3];
rz(2.5276195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7172598) q[2];
sx q[2];
rz(-0.029742664) q[2];
sx q[2];
rz(-2.7359803) q[2];
rz(2.3154216) q[3];
sx q[3];
rz(-0.043353733) q[3];
sx q[3];
rz(0.9894754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9721603) q[0];
sx q[0];
rz(-2.5763474) q[0];
sx q[0];
rz(-1.3954847) q[0];
rz(-1.5211498) q[1];
sx q[1];
rz(-2.4902159) q[1];
sx q[1];
rz(-1.7924538) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1698818) q[0];
sx q[0];
rz(-0.75877178) q[0];
sx q[0];
rz(-1.8797329) q[0];
rz(-pi) q[1];
rz(-1.3302035) q[2];
sx q[2];
rz(-1.4052011) q[2];
sx q[2];
rz(-2.9179433) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1101035) q[1];
sx q[1];
rz(-1.5564791) q[1];
sx q[1];
rz(-1.5107461) q[1];
rz(-1.5557846) q[3];
sx q[3];
rz(-1.540657) q[3];
sx q[3];
rz(2.3127401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32538432) q[2];
sx q[2];
rz(-3.1367229) q[2];
sx q[2];
rz(1.4332786) q[2];
rz(-1.2895182) q[3];
sx q[3];
rz(-3.0472445) q[3];
sx q[3];
rz(1.2963699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077944003) q[0];
sx q[0];
rz(-1.0366806) q[0];
sx q[0];
rz(0.59564577) q[0];
rz(3.0972277) q[1];
sx q[1];
rz(-2.8652419) q[1];
sx q[1];
rz(-1.2038632) q[1];
rz(-2.9246671) q[2];
sx q[2];
rz(-2.0382765) q[2];
sx q[2];
rz(-2.9983476) q[2];
rz(-0.065905215) q[3];
sx q[3];
rz(-1.1614945) q[3];
sx q[3];
rz(0.11621034) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
