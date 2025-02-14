OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.92276031) q[0];
sx q[0];
rz(-3.0781167) q[0];
sx q[0];
rz(1.3702962) q[0];
rz(0.81075794) q[1];
sx q[1];
rz(1.6914565) q[1];
sx q[1];
rz(9.20426) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0974183) q[0];
sx q[0];
rz(-0.83651354) q[0];
sx q[0];
rz(0.091139779) q[0];
x q[1];
rz(-0.96376586) q[2];
sx q[2];
rz(-2.6015601) q[2];
sx q[2];
rz(-2.2390197) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5223536) q[1];
sx q[1];
rz(-1.5818038) q[1];
sx q[1];
rz(0.22214823) q[1];
x q[2];
rz(2.9412782) q[3];
sx q[3];
rz(-1.2216) q[3];
sx q[3];
rz(-1.136245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0567131) q[2];
sx q[2];
rz(-3.1352545) q[2];
sx q[2];
rz(-0.81981266) q[2];
rz(-3.1208755) q[3];
sx q[3];
rz(-2.8262704) q[3];
sx q[3];
rz(-1.0516385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3725975) q[0];
sx q[0];
rz(-0.1854493) q[0];
sx q[0];
rz(-2.4107667) q[0];
rz(1.7164879) q[1];
sx q[1];
rz(-0.72530472) q[1];
sx q[1];
rz(-0.46754974) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0529405) q[0];
sx q[0];
rz(-3.0346617) q[0];
sx q[0];
rz(-1.8631975) q[0];
x q[1];
rz(-0.40981648) q[2];
sx q[2];
rz(-1.5404319) q[2];
sx q[2];
rz(1.3156459) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4380355) q[1];
sx q[1];
rz(-0.86895567) q[1];
sx q[1];
rz(1.7662394) q[1];
x q[2];
rz(-2.8465956) q[3];
sx q[3];
rz(-0.85912356) q[3];
sx q[3];
rz(2.0509348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1656437) q[2];
sx q[2];
rz(-0.011198137) q[2];
sx q[2];
rz(0.33640081) q[2];
rz(2.8339556) q[3];
sx q[3];
rz(-3.1308789) q[3];
sx q[3];
rz(-0.78905869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-3.0129608) q[0];
sx q[0];
rz(-2.9427981) q[0];
sx q[0];
rz(3.0096753) q[0];
rz(-1.4384653) q[1];
sx q[1];
rz(-1.7273644) q[1];
sx q[1];
rz(1.6579113) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8701829) q[0];
sx q[0];
rz(-2.462059) q[0];
sx q[0];
rz(-0.45241286) q[0];
rz(1.5533812) q[2];
sx q[2];
rz(-1.546374) q[2];
sx q[2];
rz(-1.2877854) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7805788) q[1];
sx q[1];
rz(-0.11615651) q[1];
sx q[1];
rz(-1.1067944) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0555757) q[3];
sx q[3];
rz(-1.6413851) q[3];
sx q[3];
rz(0.020078192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.57054532) q[2];
sx q[2];
rz(-1.3344301) q[2];
sx q[2];
rz(-1.6837616) q[2];
rz(-2.3633862) q[3];
sx q[3];
rz(-3.1359378) q[3];
sx q[3];
rz(0.2125423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1218162) q[0];
sx q[0];
rz(-1.3874522) q[0];
sx q[0];
rz(2.8507267) q[0];
rz(1.5930755) q[1];
sx q[1];
rz(-0.80268186) q[1];
sx q[1];
rz(-0.026550857) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0565467) q[0];
sx q[0];
rz(-1.9246058) q[0];
sx q[0];
rz(-2.7744571) q[0];
rz(-0.45088972) q[2];
sx q[2];
rz(-1.6188545) q[2];
sx q[2];
rz(-0.61651905) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8513048) q[1];
sx q[1];
rz(-1.0920224) q[1];
sx q[1];
rz(-0.7022052) q[1];
rz(-0.0031626458) q[3];
sx q[3];
rz(-1.589865) q[3];
sx q[3];
rz(0.17584383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2619541) q[2];
sx q[2];
rz(-3.1231572) q[2];
sx q[2];
rz(-0.43799841) q[2];
rz(-1.1079463) q[3];
sx q[3];
rz(-0.017939311) q[3];
sx q[3];
rz(1.3560791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22439013) q[0];
sx q[0];
rz(-0.49773911) q[0];
sx q[0];
rz(-0.14601953) q[0];
rz(-0.13305013) q[1];
sx q[1];
rz(-1.0597884) q[1];
sx q[1];
rz(0.45566794) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5815703) q[0];
sx q[0];
rz(-3.0455493) q[0];
sx q[0];
rz(-2.9314629) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5495569) q[2];
sx q[2];
rz(-1.8285654) q[2];
sx q[2];
rz(1.7532312) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4670438) q[1];
sx q[1];
rz(-1.4321126) q[1];
sx q[1];
rz(1.8013181) q[1];
rz(-pi) q[2];
rz(-1.5553873) q[3];
sx q[3];
rz(-1.3432627) q[3];
sx q[3];
rz(-1.4206545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8749775) q[2];
sx q[2];
rz(-0.019860331) q[2];
sx q[2];
rz(1.88545) q[2];
rz(-2.6202294) q[3];
sx q[3];
rz(-0.25786906) q[3];
sx q[3];
rz(2.7969587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1954023) q[0];
sx q[0];
rz(-0.38637105) q[0];
sx q[0];
rz(2.572701) q[0];
rz(-2.9733114) q[1];
sx q[1];
rz(-1.5852837) q[1];
sx q[1];
rz(3.1313484) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5880347) q[0];
sx q[0];
rz(-1.7419683) q[0];
sx q[0];
rz(1.5215015) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.019316761) q[2];
sx q[2];
rz(-1.3903119) q[2];
sx q[2];
rz(1.3767124) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.467874) q[1];
sx q[1];
rz(-1.5281902) q[1];
sx q[1];
rz(0.03454476) q[1];
x q[2];
rz(0.63416173) q[3];
sx q[3];
rz(-1.595257) q[3];
sx q[3];
rz(2.0195877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9804618) q[2];
sx q[2];
rz(-2.9176596) q[2];
sx q[2];
rz(-2.0437608) q[2];
rz(-2.791642) q[3];
sx q[3];
rz(-1.2772468) q[3];
sx q[3];
rz(0.85489517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
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
rz(-1.8304652) q[1];
sx q[1];
rz(-2.6268112) q[1];
sx q[1];
rz(0.11022551) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.742934) q[0];
sx q[0];
rz(-0.1029108) q[0];
sx q[0];
rz(-2.8714931) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1378205) q[2];
sx q[2];
rz(-1.5761778) q[2];
sx q[2];
rz(2.7710235) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4264723) q[1];
sx q[1];
rz(-2.7139152) q[1];
sx q[1];
rz(-1.5918455) q[1];
rz(-0.73506395) q[3];
sx q[3];
rz(-2.1218315) q[3];
sx q[3];
rz(-2.6363274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8608287) q[2];
sx q[2];
rz(-3.132756) q[2];
sx q[2];
rz(2.3542118) q[2];
rz(-0.27960882) q[3];
sx q[3];
rz(-0.044737261) q[3];
sx q[3];
rz(0.71981049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0019919458) q[0];
sx q[0];
rz(-2.3878492) q[0];
sx q[0];
rz(1.4308223) q[0];
rz(0.17499533) q[1];
sx q[1];
rz(-2.4038959) q[1];
sx q[1];
rz(-1.7134679) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7044107) q[0];
sx q[0];
rz(-1.5488273) q[0];
sx q[0];
rz(-0.0030513838) q[0];
rz(-pi) q[1];
rz(-1.5634086) q[2];
sx q[2];
rz(-1.566717) q[2];
sx q[2];
rz(2.5830327) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5918846) q[1];
sx q[1];
rz(-1.8401405) q[1];
sx q[1];
rz(2.9742105) q[1];
rz(-pi) q[2];
rz(2.4739856) q[3];
sx q[3];
rz(-0.35304754) q[3];
sx q[3];
rz(-1.8065344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3565107) q[2];
sx q[2];
rz(-0.76866895) q[2];
sx q[2];
rz(1.53995) q[2];
rz(-2.8421863) q[3];
sx q[3];
rz(-1.4842024) q[3];
sx q[3];
rz(-0.33777344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47141075) q[0];
sx q[0];
rz(-3.1306559) q[0];
sx q[0];
rz(0.48728824) q[0];
rz(-1.7022853) q[1];
sx q[1];
rz(-2.3479562) q[1];
sx q[1];
rz(2.7557441) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7283994) q[0];
sx q[0];
rz(-1.5110713) q[0];
sx q[0];
rz(-2.7483536) q[0];
rz(1.3273238) q[2];
sx q[2];
rz(-0.43691942) q[2];
sx q[2];
rz(-1.5658558) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.072695168) q[1];
sx q[1];
rz(-2.4494315) q[1];
sx q[1];
rz(-0.051874224) q[1];
rz(-pi) q[2];
rz(0.10200067) q[3];
sx q[3];
rz(-0.33683646) q[3];
sx q[3];
rz(-2.5276195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4243329) q[2];
sx q[2];
rz(-0.029742664) q[2];
sx q[2];
rz(2.7359803) q[2];
rz(-2.3154216) q[3];
sx q[3];
rz(-3.0982389) q[3];
sx q[3];
rz(0.9894754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9721603) q[0];
sx q[0];
rz(-2.5763474) q[0];
sx q[0];
rz(1.3954847) q[0];
rz(-1.5211498) q[1];
sx q[1];
rz(-2.4902159) q[1];
sx q[1];
rz(1.3491389) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5842218) q[0];
sx q[0];
rz(-2.2855846) q[0];
sx q[0];
rz(-2.8609311) q[0];
rz(-pi) q[1];
rz(-2.1824532) q[2];
sx q[2];
rz(-0.29116071) q[2];
sx q[2];
rz(-2.3860402) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1101035) q[1];
sx q[1];
rz(-1.5851136) q[1];
sx q[1];
rz(-1.5107461) q[1];
rz(-pi) q[2];
rz(-1.5557846) q[3];
sx q[3];
rz(-1.6009357) q[3];
sx q[3];
rz(0.82885259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8162083) q[2];
sx q[2];
rz(-3.1367229) q[2];
sx q[2];
rz(1.4332786) q[2];
rz(1.2895182) q[3];
sx q[3];
rz(-0.09434814) q[3];
sx q[3];
rz(1.2963699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-1.9738214) q[2];
sx q[2];
rz(-2.6296305) q[2];
sx q[2];
rz(0.59811022) q[2];
rz(-1.4201319) q[3];
sx q[3];
rz(-0.41427783) q[3];
sx q[3];
rz(-2.8610341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
