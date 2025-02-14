OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27999347) q[0];
sx q[0];
rz(-1.0323098) q[0];
sx q[0];
rz(1.9814459) q[0];
rz(-0.87814826) q[1];
sx q[1];
rz(-2.703009) q[1];
sx q[1];
rz(2.261472) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47326776) q[0];
sx q[0];
rz(-0.62444326) q[0];
sx q[0];
rz(-1.8683598) q[0];
x q[1];
rz(-1.7574667) q[2];
sx q[2];
rz(-2.7323033) q[2];
sx q[2];
rz(-3.0011506) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9701201) q[1];
sx q[1];
rz(-0.62343854) q[1];
sx q[1];
rz(1.2894676) q[1];
rz(-pi) q[2];
rz(1.6900469) q[3];
sx q[3];
rz(-1.1195983) q[3];
sx q[3];
rz(0.50332068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0264414) q[2];
sx q[2];
rz(-0.5618962) q[2];
sx q[2];
rz(-1.119841) q[2];
rz(-2.1378873) q[3];
sx q[3];
rz(-0.34990889) q[3];
sx q[3];
rz(0.85163918) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.411946) q[0];
sx q[0];
rz(-0.8929407) q[0];
sx q[0];
rz(-0.44806421) q[0];
rz(1.0753151) q[1];
sx q[1];
rz(-1.0766462) q[1];
sx q[1];
rz(-3.101128) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4141636) q[0];
sx q[0];
rz(-2.232904) q[0];
sx q[0];
rz(0.24358769) q[0];
x q[1];
rz(-1.1574817) q[2];
sx q[2];
rz(-1.6530452) q[2];
sx q[2];
rz(0.27378191) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.65461797) q[1];
sx q[1];
rz(-0.27187133) q[1];
sx q[1];
rz(-0.15403037) q[1];
rz(0.31261698) q[3];
sx q[3];
rz(-1.5789869) q[3];
sx q[3];
rz(2.9401478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2761817) q[2];
sx q[2];
rz(-2.2409596) q[2];
sx q[2];
rz(2.6382228) q[2];
rz(1.6940176) q[3];
sx q[3];
rz(-1.9083128) q[3];
sx q[3];
rz(-2.2342822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.68509787) q[0];
sx q[0];
rz(-0.95360294) q[0];
sx q[0];
rz(2.8431235) q[0];
rz(-1.5553156) q[1];
sx q[1];
rz(-1.7319873) q[1];
sx q[1];
rz(1.635199) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28373805) q[0];
sx q[0];
rz(-1.8687792) q[0];
sx q[0];
rz(-1.6493173) q[0];
rz(-2.2573264) q[2];
sx q[2];
rz(-1.7873014) q[2];
sx q[2];
rz(-0.51196438) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.071722592) q[1];
sx q[1];
rz(-0.89080526) q[1];
sx q[1];
rz(-1.7417439) q[1];
rz(-pi) q[2];
rz(-1.7603586) q[3];
sx q[3];
rz(-1.624057) q[3];
sx q[3];
rz(1.9587751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7895268) q[2];
sx q[2];
rz(-2.3730998) q[2];
sx q[2];
rz(0.7716158) q[2];
rz(1.4113034) q[3];
sx q[3];
rz(-1.227042) q[3];
sx q[3];
rz(-0.82967776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70374933) q[0];
sx q[0];
rz(-2.4898536) q[0];
sx q[0];
rz(1.334345) q[0];
rz(1.5127381) q[1];
sx q[1];
rz(-1.6213497) q[1];
sx q[1];
rz(0.13660647) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.893394) q[0];
sx q[0];
rz(-1.8780439) q[0];
sx q[0];
rz(0.91715468) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0380546) q[2];
sx q[2];
rz(-1.8788726) q[2];
sx q[2];
rz(-3.0335226) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1699511) q[1];
sx q[1];
rz(-1.882269) q[1];
sx q[1];
rz(0.71572742) q[1];
rz(-0.81035536) q[3];
sx q[3];
rz(-2.6763066) q[3];
sx q[3];
rz(-0.76818675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2887349) q[2];
sx q[2];
rz(-1.4790269) q[2];
sx q[2];
rz(0.23957254) q[2];
rz(-2.1824956) q[3];
sx q[3];
rz(-1.164091) q[3];
sx q[3];
rz(2.1036928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94944823) q[0];
sx q[0];
rz(-2.0862155) q[0];
sx q[0];
rz(1.5437833) q[0];
rz(-1.1100618) q[1];
sx q[1];
rz(-1.9381356) q[1];
sx q[1];
rz(1.4470626) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8657762) q[0];
sx q[0];
rz(-1.5636684) q[0];
sx q[0];
rz(2.843186) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1008312) q[2];
sx q[2];
rz(-1.9136179) q[2];
sx q[2];
rz(2.0343604) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.93433468) q[1];
sx q[1];
rz(-0.90597766) q[1];
sx q[1];
rz(-1.7938754) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5041833) q[3];
sx q[3];
rz(-0.64475179) q[3];
sx q[3];
rz(-1.7774323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1453104) q[2];
sx q[2];
rz(-1.1386394) q[2];
sx q[2];
rz(0.81926695) q[2];
rz(2.1938358) q[3];
sx q[3];
rz(-2.3534333) q[3];
sx q[3];
rz(1.88571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6325697) q[0];
sx q[0];
rz(-3.0036354) q[0];
sx q[0];
rz(0.55009681) q[0];
rz(-2.8944648) q[1];
sx q[1];
rz(-1.1526266) q[1];
sx q[1];
rz(0.94608847) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6912963) q[0];
sx q[0];
rz(-1.2142203) q[0];
sx q[0];
rz(-2.8647115) q[0];
rz(1.7839512) q[2];
sx q[2];
rz(-1.5098796) q[2];
sx q[2];
rz(-2.3062381) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.94080776) q[1];
sx q[1];
rz(-2.0403123) q[1];
sx q[1];
rz(1.5307309) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9027418) q[3];
sx q[3];
rz(-1.4689886) q[3];
sx q[3];
rz(-0.073697173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.75521022) q[2];
sx q[2];
rz(-2.5911665) q[2];
sx q[2];
rz(2.4605301) q[2];
rz(-0.81124535) q[3];
sx q[3];
rz(-1.0439876) q[3];
sx q[3];
rz(0.55027858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1193467) q[0];
sx q[0];
rz(-2.5211054) q[0];
sx q[0];
rz(-2.3959809) q[0];
rz(1.9781808) q[1];
sx q[1];
rz(-2.5201576) q[1];
sx q[1];
rz(2.9643639) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39800571) q[0];
sx q[0];
rz(-1.9015802) q[0];
sx q[0];
rz(2.8547006) q[0];
x q[1];
rz(1.4800328) q[2];
sx q[2];
rz(-1.5077733) q[2];
sx q[2];
rz(1.3985967) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.16645277) q[1];
sx q[1];
rz(-1.9362294) q[1];
sx q[1];
rz(2.1154984) q[1];
rz(2.4468462) q[3];
sx q[3];
rz(-1.5383895) q[3];
sx q[3];
rz(2.29914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.54997286) q[2];
sx q[2];
rz(-2.1592996) q[2];
sx q[2];
rz(-1.4493235) q[2];
rz(-0.6226522) q[3];
sx q[3];
rz(-2.6144274) q[3];
sx q[3];
rz(-1.8221633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016668884) q[0];
sx q[0];
rz(-1.387384) q[0];
sx q[0];
rz(1.6296052) q[0];
rz(-2.2177057) q[1];
sx q[1];
rz(-1.7216564) q[1];
sx q[1];
rz(0.61663827) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9074207) q[0];
sx q[0];
rz(-1.7230095) q[0];
sx q[0];
rz(-0.9665177) q[0];
rz(-pi) q[1];
rz(-2.6737164) q[2];
sx q[2];
rz(-1.6098002) q[2];
sx q[2];
rz(2.6371961) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.42378572) q[1];
sx q[1];
rz(-0.95342481) q[1];
sx q[1];
rz(-1.6758534) q[1];
x q[2];
rz(-0.40469669) q[3];
sx q[3];
rz(-2.5202978) q[3];
sx q[3];
rz(2.9018847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7743249) q[2];
sx q[2];
rz(-0.35949817) q[2];
sx q[2];
rz(-0.16505879) q[2];
rz(2.3955691) q[3];
sx q[3];
rz(-1.5187998) q[3];
sx q[3];
rz(0.042796854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5999254) q[0];
sx q[0];
rz(-0.2921108) q[0];
sx q[0];
rz(0.41411972) q[0];
rz(-2.7060624) q[1];
sx q[1];
rz(-2.6256517) q[1];
sx q[1];
rz(-1.863265) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3992452) q[0];
sx q[0];
rz(-1.81184) q[0];
sx q[0];
rz(2.9910415) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8299471) q[2];
sx q[2];
rz(-1.4053182) q[2];
sx q[2];
rz(-2.8143425) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8202701) q[1];
sx q[1];
rz(-1.6386809) q[1];
sx q[1];
rz(1.5826788) q[1];
rz(-pi) q[2];
rz(-2.0986365) q[3];
sx q[3];
rz(-1.5164638) q[3];
sx q[3];
rz(-2.899037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4182959) q[2];
sx q[2];
rz(-1.7593242) q[2];
sx q[2];
rz(2.3988775) q[2];
rz(-2.3501979) q[3];
sx q[3];
rz(-2.4580038) q[3];
sx q[3];
rz(1.7323823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0940229) q[0];
sx q[0];
rz(-0.3293193) q[0];
sx q[0];
rz(-0.91162115) q[0];
rz(2.1564663) q[1];
sx q[1];
rz(-0.42675012) q[1];
sx q[1];
rz(1.8034579) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0165923) q[0];
sx q[0];
rz(-0.58412281) q[0];
sx q[0];
rz(0.97726269) q[0];
x q[1];
rz(0.31735898) q[2];
sx q[2];
rz(-1.9322763) q[2];
sx q[2];
rz(-1.1731479) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.12305189) q[1];
sx q[1];
rz(-1.9705074) q[1];
sx q[1];
rz(-2.2247061) q[1];
rz(-pi) q[2];
rz(1.1971522) q[3];
sx q[3];
rz(-1.8836397) q[3];
sx q[3];
rz(-2.2598303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7868339) q[2];
sx q[2];
rz(-1.6666731) q[2];
sx q[2];
rz(-1.9762529) q[2];
rz(1.3960086) q[3];
sx q[3];
rz(-1.952012) q[3];
sx q[3];
rz(1.5918119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6373445) q[0];
sx q[0];
rz(-1.5865542) q[0];
sx q[0];
rz(-2.4028548) q[0];
rz(-2.5042116) q[1];
sx q[1];
rz(-1.6045872) q[1];
sx q[1];
rz(2.3789023) q[1];
rz(0.1461045) q[2];
sx q[2];
rz(-1.1510104) q[2];
sx q[2];
rz(0.83465464) q[2];
rz(2.4828667) q[3];
sx q[3];
rz(-0.6686419) q[3];
sx q[3];
rz(1.1262459) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
