OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.926446) q[0];
sx q[0];
rz(-0.67987052) q[0];
sx q[0];
rz(3.0409066) q[0];
rz(-2.6868532) q[1];
sx q[1];
rz(-2.4603031) q[1];
sx q[1];
rz(1.0256306) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68368841) q[0];
sx q[0];
rz(-1.6207028) q[0];
sx q[0];
rz(-1.7662952) q[0];
x q[1];
rz(-1.6221912) q[2];
sx q[2];
rz(-1.8784461) q[2];
sx q[2];
rz(-0.6531229) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6903578) q[1];
sx q[1];
rz(-0.75690311) q[1];
sx q[1];
rz(-3.0249658) q[1];
x q[2];
rz(1.0354393) q[3];
sx q[3];
rz(-0.97781721) q[3];
sx q[3];
rz(1.6361332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.92496282) q[2];
sx q[2];
rz(-1.7731885) q[2];
sx q[2];
rz(1.487757) q[2];
rz(1.9803068) q[3];
sx q[3];
rz(-1.0043283) q[3];
sx q[3];
rz(-2.0465093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0682812) q[0];
sx q[0];
rz(-0.40638766) q[0];
sx q[0];
rz(0.43346369) q[0];
rz(-2.0794226) q[1];
sx q[1];
rz(-1.559161) q[1];
sx q[1];
rz(2.4210222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5262289) q[0];
sx q[0];
rz(-3.0096292) q[0];
sx q[0];
rz(1.855164) q[0];
rz(-pi) q[1];
rz(2.8938204) q[2];
sx q[2];
rz(-2.0907913) q[2];
sx q[2];
rz(1.2733851) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7097672) q[1];
sx q[1];
rz(-2.7569237) q[1];
sx q[1];
rz(-2.7363214) q[1];
rz(-pi) q[2];
rz(2.0907164) q[3];
sx q[3];
rz(-2.5934016) q[3];
sx q[3];
rz(-0.91778008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.34708193) q[2];
sx q[2];
rz(-1.2257267) q[2];
sx q[2];
rz(-2.1573055) q[2];
rz(2.2847564) q[3];
sx q[3];
rz(-2.555116) q[3];
sx q[3];
rz(1.2543359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03265753) q[0];
sx q[0];
rz(-0.28626838) q[0];
sx q[0];
rz(2.3987067) q[0];
rz(-0.0097097857) q[1];
sx q[1];
rz(-1.5747036) q[1];
sx q[1];
rz(0.28365338) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24697248) q[0];
sx q[0];
rz(-0.63480154) q[0];
sx q[0];
rz(-0.60946861) q[0];
x q[1];
rz(2.4667612) q[2];
sx q[2];
rz(-0.25383224) q[2];
sx q[2];
rz(0.93364894) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4176826) q[1];
sx q[1];
rz(-1.4028478) q[1];
sx q[1];
rz(-1.8803094) q[1];
rz(-0.023223485) q[3];
sx q[3];
rz(-1.4693714) q[3];
sx q[3];
rz(-1.9254799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0297086) q[2];
sx q[2];
rz(-1.6969029) q[2];
sx q[2];
rz(-0.65579826) q[2];
rz(0.2119952) q[3];
sx q[3];
rz(-0.63181221) q[3];
sx q[3];
rz(1.1727715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1126605) q[0];
sx q[0];
rz(-2.6472968) q[0];
sx q[0];
rz(-1.5546881) q[0];
rz(-0.2298062) q[1];
sx q[1];
rz(-0.72139144) q[1];
sx q[1];
rz(-1.302964) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.479871) q[0];
sx q[0];
rz(-2.4085143) q[0];
sx q[0];
rz(0.69723155) q[0];
rz(-pi) q[1];
rz(2.8726758) q[2];
sx q[2];
rz(-1.285621) q[2];
sx q[2];
rz(1.1095003) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.069217056) q[1];
sx q[1];
rz(-1.0472968) q[1];
sx q[1];
rz(1.5956466) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0082208) q[3];
sx q[3];
rz(-1.2689212) q[3];
sx q[3];
rz(0.20193298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.135123) q[2];
sx q[2];
rz(-2.3389356) q[2];
sx q[2];
rz(-2.6378677) q[2];
rz(-1.6040365) q[3];
sx q[3];
rz(-1.1732912) q[3];
sx q[3];
rz(1.5776207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2876494) q[0];
sx q[0];
rz(-0.33677736) q[0];
sx q[0];
rz(-0.42959282) q[0];
rz(2.8751539) q[1];
sx q[1];
rz(-2.3315505) q[1];
sx q[1];
rz(-1.0688759) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3392068) q[0];
sx q[0];
rz(-0.9046208) q[0];
sx q[0];
rz(-2.1016913) q[0];
rz(-pi) q[1];
rz(-0.30854723) q[2];
sx q[2];
rz(-0.82627901) q[2];
sx q[2];
rz(-2.770785) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4539448) q[1];
sx q[1];
rz(-1.8275211) q[1];
sx q[1];
rz(2.686108) q[1];
x q[2];
rz(-0.74917082) q[3];
sx q[3];
rz(-0.35532829) q[3];
sx q[3];
rz(-0.65566777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8916919) q[2];
sx q[2];
rz(-2.4155858) q[2];
sx q[2];
rz(-2.967584) q[2];
rz(1.8788667) q[3];
sx q[3];
rz(-1.5401253) q[3];
sx q[3];
rz(-1.5549829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0083017666) q[0];
sx q[0];
rz(-2.1228078) q[0];
sx q[0];
rz(0.89023501) q[0];
rz(1.7091735) q[1];
sx q[1];
rz(-2.1232257) q[1];
sx q[1];
rz(0.11996809) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9667119) q[0];
sx q[0];
rz(-0.57583664) q[0];
sx q[0];
rz(2.8737646) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2416297) q[2];
sx q[2];
rz(-0.53027486) q[2];
sx q[2];
rz(1.2880304) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6585582) q[1];
sx q[1];
rz(-2.2081828) q[1];
sx q[1];
rz(-1.263522) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31474416) q[3];
sx q[3];
rz(-0.56296825) q[3];
sx q[3];
rz(-1.1551577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1098108) q[2];
sx q[2];
rz(-1.4795156) q[2];
sx q[2];
rz(-0.39153448) q[2];
rz(1.7560962) q[3];
sx q[3];
rz(-2.9031495) q[3];
sx q[3];
rz(2.5029206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2077797) q[0];
sx q[0];
rz(-0.55140984) q[0];
sx q[0];
rz(-1.3054003) q[0];
rz(0.50453672) q[1];
sx q[1];
rz(-1.4089818) q[1];
sx q[1];
rz(-2.9171468) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2674057) q[0];
sx q[0];
rz(-1.5541557) q[0];
sx q[0];
rz(-1.5865302) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7573757) q[2];
sx q[2];
rz(-2.6165892) q[2];
sx q[2];
rz(0.25431654) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9516884) q[1];
sx q[1];
rz(-2.7153141) q[1];
sx q[1];
rz(-0.87949068) q[1];
rz(-pi) q[2];
rz(2.283758) q[3];
sx q[3];
rz(-2.668475) q[3];
sx q[3];
rz(-0.83698646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.18181431) q[2];
sx q[2];
rz(-1.10428) q[2];
sx q[2];
rz(-2.1941198) q[2];
rz(0.13253658) q[3];
sx q[3];
rz(-2.4819076) q[3];
sx q[3];
rz(-2.1555677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.9082311) q[0];
sx q[0];
rz(-2.7955604) q[0];
sx q[0];
rz(0.37286266) q[0];
rz(2.8064959) q[1];
sx q[1];
rz(-1.9357888) q[1];
sx q[1];
rz(-2.0705409) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.459254) q[0];
sx q[0];
rz(-1.5821348) q[0];
sx q[0];
rz(0.069123231) q[0];
rz(-0.33800563) q[2];
sx q[2];
rz(-2.5635984) q[2];
sx q[2];
rz(-1.6737311) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.11406273) q[1];
sx q[1];
rz(-2.8004871) q[1];
sx q[1];
rz(-2.3807664) q[1];
rz(-pi) q[2];
rz(-1.8926431) q[3];
sx q[3];
rz(-1.6610356) q[3];
sx q[3];
rz(-2.7863996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7143453) q[2];
sx q[2];
rz(-1.5421966) q[2];
sx q[2];
rz(2.074312) q[2];
rz(-0.37211564) q[3];
sx q[3];
rz(-0.99388638) q[3];
sx q[3];
rz(-0.20017643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.9099971) q[0];
sx q[0];
rz(-0.59760004) q[0];
sx q[0];
rz(-1.4071314) q[0];
rz(-1.0435957) q[1];
sx q[1];
rz(-2.2035972) q[1];
sx q[1];
rz(-2.6503906) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3160014) q[0];
sx q[0];
rz(-1.2818616) q[0];
sx q[0];
rz(-0.87744539) q[0];
x q[1];
rz(-1.4004854) q[2];
sx q[2];
rz(-1.7372157) q[2];
sx q[2];
rz(0.8569878) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.72043967) q[1];
sx q[1];
rz(-0.3657839) q[1];
sx q[1];
rz(0.37558742) q[1];
rz(-0.37685412) q[3];
sx q[3];
rz(-1.4026733) q[3];
sx q[3];
rz(-0.73790077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4444943) q[2];
sx q[2];
rz(-1.5105057) q[2];
sx q[2];
rz(2.0972882) q[2];
rz(0.77110428) q[3];
sx q[3];
rz(-1.5325357) q[3];
sx q[3];
rz(-0.3573629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.20030178) q[0];
sx q[0];
rz(-1.8481978) q[0];
sx q[0];
rz(-1.6218761) q[0];
rz(1.3007851) q[1];
sx q[1];
rz(-1.8404574) q[1];
sx q[1];
rz(2.4470952) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3169169) q[0];
sx q[0];
rz(-1.3429214) q[0];
sx q[0];
rz(-3.1275463) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7037665) q[2];
sx q[2];
rz(-2.0914234) q[2];
sx q[2];
rz(0.43308738) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2934781) q[1];
sx q[1];
rz(-1.2367931) q[1];
sx q[1];
rz(-0.011718463) q[1];
rz(-2.3894074) q[3];
sx q[3];
rz(-1.859741) q[3];
sx q[3];
rz(2.1352701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0684315) q[2];
sx q[2];
rz(-1.7794926) q[2];
sx q[2];
rz(0.51363242) q[2];
rz(-2.8237776) q[3];
sx q[3];
rz(-1.603926) q[3];
sx q[3];
rz(2.1224799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1401055) q[0];
sx q[0];
rz(-1.6248063) q[0];
sx q[0];
rz(-0.90909062) q[0];
rz(-1.2263251) q[1];
sx q[1];
rz(-1.6188123) q[1];
sx q[1];
rz(0.34368044) q[1];
rz(0.60811715) q[2];
sx q[2];
rz(-0.79610745) q[2];
sx q[2];
rz(2.983762) q[2];
rz(-1.0448643) q[3];
sx q[3];
rz(-1.5689701) q[3];
sx q[3];
rz(0.29472385) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
