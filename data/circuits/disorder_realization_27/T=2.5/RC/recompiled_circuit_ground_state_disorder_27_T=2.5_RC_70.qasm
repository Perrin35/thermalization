OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5840924) q[0];
sx q[0];
rz(3.1182365) q[0];
sx q[0];
rz(11.630848) q[0];
rz(1.525653) q[1];
sx q[1];
rz(4.6620044) q[1];
sx q[1];
rz(9.6992156) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0204791) q[0];
sx q[0];
rz(-1.6389567) q[0];
sx q[0];
rz(-1.6421374) q[0];
rz(-pi) q[1];
rz(2.0361314) q[2];
sx q[2];
rz(-2.0623775) q[2];
sx q[2];
rz(0.21869379) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5554065) q[1];
sx q[1];
rz(-1.6073391) q[1];
sx q[1];
rz(-0.010374109) q[1];
rz(-pi) q[2];
rz(-0.07368659) q[3];
sx q[3];
rz(-1.4006124) q[3];
sx q[3];
rz(-0.88995349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9258257) q[2];
sx q[2];
rz(-0.01130686) q[2];
sx q[2];
rz(-2.0634148) q[2];
rz(2.3128541) q[3];
sx q[3];
rz(-1.6308035) q[3];
sx q[3];
rz(-0.80882788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082212903) q[0];
sx q[0];
rz(-1.2780715) q[0];
sx q[0];
rz(1.7510121) q[0];
rz(-1.4311721) q[1];
sx q[1];
rz(-3.1371959) q[1];
sx q[1];
rz(1.4336525) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2947049) q[0];
sx q[0];
rz(-1.6754912) q[0];
sx q[0];
rz(-2.2082735) q[0];
rz(-pi) q[1];
rz(1.6054691) q[2];
sx q[2];
rz(-0.44402396) q[2];
sx q[2];
rz(1.6426067) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8535276) q[1];
sx q[1];
rz(-1.57987) q[1];
sx q[1];
rz(0.00041043343) q[1];
x q[2];
rz(1.6363899) q[3];
sx q[3];
rz(-2.3573313) q[3];
sx q[3];
rz(-1.0515056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.37890515) q[2];
sx q[2];
rz(-1.5964369) q[2];
sx q[2];
rz(-1.5691266) q[2];
rz(0.76111859) q[3];
sx q[3];
rz(-3.0877536) q[3];
sx q[3];
rz(-1.9332473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91956562) q[0];
sx q[0];
rz(-0.61028218) q[0];
sx q[0];
rz(0.56184226) q[0];
rz(1.5737083) q[1];
sx q[1];
rz(-1.5627197) q[1];
sx q[1];
rz(3.124253) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8652788) q[0];
sx q[0];
rz(-1.9863578) q[0];
sx q[0];
rz(2.1091976) q[0];
rz(-pi) q[1];
rz(2.1465535) q[2];
sx q[2];
rz(-0.77313609) q[2];
sx q[2];
rz(-2.7553158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8120809) q[1];
sx q[1];
rz(-1.2083665) q[1];
sx q[1];
rz(3.1245438) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51029499) q[3];
sx q[3];
rz(-1.5908352) q[3];
sx q[3];
rz(0.35773548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5138381) q[2];
sx q[2];
rz(-1.4189812) q[2];
sx q[2];
rz(-0.55297744) q[2];
rz(-1.9364457) q[3];
sx q[3];
rz(-1.5670992) q[3];
sx q[3];
rz(-1.5945826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7503081) q[0];
sx q[0];
rz(-2.1109844) q[0];
sx q[0];
rz(-1.7872101) q[0];
rz(-1.7693819) q[1];
sx q[1];
rz(-0.0028227614) q[1];
sx q[1];
rz(-1.7599958) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1005122) q[0];
sx q[0];
rz(-1.1100475) q[0];
sx q[0];
rz(-0.45188388) q[0];
x q[1];
rz(-0.69475485) q[2];
sx q[2];
rz(-0.0036029795) q[2];
sx q[2];
rz(0.82243332) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.18455566) q[1];
sx q[1];
rz(-2.364675) q[1];
sx q[1];
rz(2.0691815) q[1];
x q[2];
rz(0.28691157) q[3];
sx q[3];
rz(-2.3182456) q[3];
sx q[3];
rz(-2.8475743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0682721) q[2];
sx q[2];
rz(-0.019651532) q[2];
sx q[2];
rz(1.2400631) q[2];
rz(-0.23010075) q[3];
sx q[3];
rz(-3.1374044) q[3];
sx q[3];
rz(-2.7219462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0536026) q[0];
sx q[0];
rz(-2.3542861) q[0];
sx q[0];
rz(1.2552274) q[0];
rz(-3.1322196) q[1];
sx q[1];
rz(-1.7724937) q[1];
sx q[1];
rz(-0.031551687) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2627336) q[0];
sx q[0];
rz(-0.53476214) q[0];
sx q[0];
rz(1.4654903) q[0];
rz(-1.2826233) q[2];
sx q[2];
rz(-2.8409578) q[2];
sx q[2];
rz(-0.86821454) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7158791) q[1];
sx q[1];
rz(-1.3495933) q[1];
sx q[1];
rz(2.1137266) q[1];
rz(-pi) q[2];
x q[2];
rz(2.191014) q[3];
sx q[3];
rz(-1.0744541) q[3];
sx q[3];
rz(-2.2994808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8072529) q[2];
sx q[2];
rz(-0.0060609239) q[2];
sx q[2];
rz(-2.3414211) q[2];
rz(2.3470894) q[3];
sx q[3];
rz(-0.032363351) q[3];
sx q[3];
rz(-0.86875027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.054258) q[0];
sx q[0];
rz(-2.9825409) q[0];
sx q[0];
rz(-1.4999088) q[0];
rz(-0.17290393) q[1];
sx q[1];
rz(-0.042363107) q[1];
sx q[1];
rz(-0.049887966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099012696) q[0];
sx q[0];
rz(-0.96262162) q[0];
sx q[0];
rz(0.30607589) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3388145) q[2];
sx q[2];
rz(-1.652431) q[2];
sx q[2];
rz(0.25914524) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8356925) q[1];
sx q[1];
rz(-0.75961194) q[1];
sx q[1];
rz(-1.2051969) q[1];
rz(-2.6088151) q[3];
sx q[3];
rz(-0.92255892) q[3];
sx q[3];
rz(-1.6870354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40338966) q[2];
sx q[2];
rz(-3.093779) q[2];
sx q[2];
rz(-1.8955463) q[2];
rz(-1.7854779) q[3];
sx q[3];
rz(-0.035364371) q[3];
sx q[3];
rz(1.7252007) q[3];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4288915) q[0];
sx q[0];
rz(-0.8242979) q[0];
sx q[0];
rz(-1.7190546) q[0];
rz(1.9046344) q[1];
sx q[1];
rz(-0.037038602) q[1];
sx q[1];
rz(2.9388156) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13064676) q[0];
sx q[0];
rz(-1.2460862) q[0];
sx q[0];
rz(2.3492011) q[0];
rz(-pi) q[1];
rz(0.76272623) q[2];
sx q[2];
rz(-0.87050754) q[2];
sx q[2];
rz(-1.3415847) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8494871) q[1];
sx q[1];
rz(-1.6260481) q[1];
sx q[1];
rz(-2.8007052) q[1];
rz(-1.0596541) q[3];
sx q[3];
rz(-1.4333925) q[3];
sx q[3];
rz(0.58611548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5286336) q[2];
sx q[2];
rz(-3.0411159) q[2];
sx q[2];
rz(2.4670777) q[2];
rz(-1.7965192) q[3];
sx q[3];
rz(-2.9967872) q[3];
sx q[3];
rz(1.323918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26789185) q[0];
sx q[0];
rz(-0.75650263) q[0];
sx q[0];
rz(0.85195136) q[0];
rz(0.1918699) q[1];
sx q[1];
rz(-3.1286897) q[1];
sx q[1];
rz(-2.875944) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11394792) q[0];
sx q[0];
rz(-2.6922604) q[0];
sx q[0];
rz(0.41037257) q[0];
rz(-pi) q[1];
rz(2.526409) q[2];
sx q[2];
rz(-1.8447723) q[2];
sx q[2];
rz(0.09466234) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3959256) q[1];
sx q[1];
rz(-0.10217459) q[1];
sx q[1];
rz(0.61421444) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1962131) q[3];
sx q[3];
rz(-1.0455971) q[3];
sx q[3];
rz(-1.7708667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77136451) q[2];
sx q[2];
rz(-3.0493272) q[2];
sx q[2];
rz(-0.51221687) q[2];
rz(-0.15906119) q[3];
sx q[3];
rz(-0.03511196) q[3];
sx q[3];
rz(-1.4353282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8856186) q[0];
sx q[0];
rz(-1.4775448) q[0];
sx q[0];
rz(1.0783827) q[0];
rz(-1.501561) q[1];
sx q[1];
rz(-2.9402132) q[1];
sx q[1];
rz(-1.5837502) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013951741) q[0];
sx q[0];
rz(-1.5553345) q[0];
sx q[0];
rz(-0.68711335) q[0];
x q[1];
rz(0.89084004) q[2];
sx q[2];
rz(-1.250923) q[2];
sx q[2];
rz(0.25242463) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.16405995) q[1];
sx q[1];
rz(-1.5674453) q[1];
sx q[1];
rz(-3.1195669) q[1];
x q[2];
rz(1.1224062) q[3];
sx q[3];
rz(-2.4351154) q[3];
sx q[3];
rz(-2.024533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.25643361) q[2];
sx q[2];
rz(-0.014545518) q[2];
sx q[2];
rz(3.0719768) q[2];
rz(0.066667892) q[3];
sx q[3];
rz(-1.0144517) q[3];
sx q[3];
rz(-2.4623509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.2081864) q[0];
sx q[0];
rz(-1.2766159) q[0];
sx q[0];
rz(0.25190121) q[0];
rz(-1.6690286) q[1];
sx q[1];
rz(-0.21569574) q[1];
sx q[1];
rz(3.0676945) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2011482) q[0];
sx q[0];
rz(-1.6415039) q[0];
sx q[0];
rz(1.9364249) q[0];
rz(-0.017357512) q[2];
sx q[2];
rz(-1.6021172) q[2];
sx q[2];
rz(0.90085122) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.628129) q[1];
sx q[1];
rz(-0.76490116) q[1];
sx q[1];
rz(2.7897863) q[1];
rz(0.68040397) q[3];
sx q[3];
rz(-2.8814253) q[3];
sx q[3];
rz(-2.5870067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4613688) q[2];
sx q[2];
rz(-3.1344487) q[2];
sx q[2];
rz(-0.76772493) q[2];
rz(1.7420306) q[3];
sx q[3];
rz(-3.1407686) q[3];
sx q[3];
rz(-0.50423938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5117699) q[0];
sx q[0];
rz(-0.98291021) q[0];
sx q[0];
rz(1.7194189) q[0];
rz(-0.024367532) q[1];
sx q[1];
rz(-2.9822646) q[1];
sx q[1];
rz(0.23039625) q[1];
rz(0.89543912) q[2];
sx q[2];
rz(-1.2661305) q[2];
sx q[2];
rz(0.2998395) q[2];
rz(-0.21655269) q[3];
sx q[3];
rz(-0.42600507) q[3];
sx q[3];
rz(1.101839) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
