OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9159311) q[0];
sx q[0];
rz(-0.8684648) q[0];
sx q[0];
rz(2.948569) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(3.5715754) q[1];
sx q[1];
rz(6.9663098) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8773168) q[0];
sx q[0];
rz(-1.5125934) q[0];
sx q[0];
rz(0.067406128) q[0];
x q[1];
rz(-1.921466) q[2];
sx q[2];
rz(-1.3436683) q[2];
sx q[2];
rz(-0.35958689) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7072304) q[1];
sx q[1];
rz(-1.1303567) q[1];
sx q[1];
rz(-0.67727725) q[1];
x q[2];
rz(1.0115511) q[3];
sx q[3];
rz(-0.64239255) q[3];
sx q[3];
rz(2.0410283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1258939) q[2];
sx q[2];
rz(-1.761972) q[2];
sx q[2];
rz(-1.0985628) q[2];
rz(2.0627608) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(-2.0400955) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1612448) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(-0.20794491) q[0];
rz(2.5646599) q[1];
sx q[1];
rz(-0.88795841) q[1];
sx q[1];
rz(-1.6764486) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10509051) q[0];
sx q[0];
rz(-2.4173792) q[0];
sx q[0];
rz(3.1378531) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25603489) q[2];
sx q[2];
rz(-1.0962152) q[2];
sx q[2];
rz(1.408996) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.84491731) q[1];
sx q[1];
rz(-2.1293318) q[1];
sx q[1];
rz(1.0122453) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74554262) q[3];
sx q[3];
rz(-1.3127483) q[3];
sx q[3];
rz(0.52449709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1229822) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(3.0318276) q[2];
rz(-0.62260735) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(-1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96317545) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(2.6254568) q[0];
rz(-2.5667045) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(2.3410472) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3886867) q[0];
sx q[0];
rz(-1.9770925) q[0];
sx q[0];
rz(-2.9857062) q[0];
rz(-pi) q[1];
rz(0.86513743) q[2];
sx q[2];
rz(-1.0312928) q[2];
sx q[2];
rz(-1.7973961) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2785981) q[1];
sx q[1];
rz(-1.9837712) q[1];
sx q[1];
rz(-1.5652656) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9242026) q[3];
sx q[3];
rz(-2.4089775) q[3];
sx q[3];
rz(-2.4917847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67733726) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(-1.7819972) q[2];
rz(0.96757403) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(-1.7165855) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1490705) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(-0.46491369) q[0];
rz(-2.7930296) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(-2.0565313) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0264498) q[0];
sx q[0];
rz(-2.0949674) q[0];
sx q[0];
rz(-1.8379704) q[0];
rz(-2.1539139) q[2];
sx q[2];
rz(-2.7104125) q[2];
sx q[2];
rz(-2.0476066) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3241868) q[1];
sx q[1];
rz(-1.4293912) q[1];
sx q[1];
rz(1.4694587) q[1];
rz(-0.61120175) q[3];
sx q[3];
rz(-1.9107606) q[3];
sx q[3];
rz(-1.1838278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.00502914) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(-2.4604649) q[2];
rz(-0.51182169) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(-0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27914771) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(1.7329247) q[0];
rz(-2.7092343) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(0.98168215) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78445804) q[0];
sx q[0];
rz(-2.031209) q[0];
sx q[0];
rz(-2.5255192) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9972516) q[2];
sx q[2];
rz(-2.855636) q[2];
sx q[2];
rz(2.2821102) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2770734) q[1];
sx q[1];
rz(-0.48950567) q[1];
sx q[1];
rz(2.8237052) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97335191) q[3];
sx q[3];
rz(-1.5095599) q[3];
sx q[3];
rz(-2.6254568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0354054) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(0.48941082) q[2];
rz(2.126746) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(1.813252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6569825) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(-0.37242517) q[0];
rz(-1.3308446) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(2.9763124) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3344791) q[0];
sx q[0];
rz(-1.4757336) q[0];
sx q[0];
rz(1.6001742) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68768244) q[2];
sx q[2];
rz(-1.624794) q[2];
sx q[2];
rz(2.8962367) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.936441) q[1];
sx q[1];
rz(-1.1679822) q[1];
sx q[1];
rz(0.58516296) q[1];
rz(-1.1092471) q[3];
sx q[3];
rz(-2.647532) q[3];
sx q[3];
rz(2.9500614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2281987) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(-1.6983263) q[2];
rz(2.7741487) q[3];
sx q[3];
rz(-1.2747217) q[3];
sx q[3];
rz(0.2120367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3570324) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(0.34926397) q[0];
rz(0.7473942) q[1];
sx q[1];
rz(-0.29577297) q[1];
sx q[1];
rz(0.73648891) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.526028) q[0];
sx q[0];
rz(-3.0897339) q[0];
sx q[0];
rz(1.2349013) q[0];
rz(-1.9095124) q[2];
sx q[2];
rz(-1.806353) q[2];
sx q[2];
rz(-0.92600694) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6474825) q[1];
sx q[1];
rz(-0.78299114) q[1];
sx q[1];
rz(-0.25581911) q[1];
x q[2];
rz(2.6408259) q[3];
sx q[3];
rz(-1.9541249) q[3];
sx q[3];
rz(-2.4367743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0050469) q[2];
sx q[2];
rz(-1.2282635) q[2];
sx q[2];
rz(2.6100256) q[2];
rz(0.68938869) q[3];
sx q[3];
rz(-1.6770984) q[3];
sx q[3];
rz(-1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625967) q[0];
sx q[0];
rz(-1.2051219) q[0];
sx q[0];
rz(1.8435562) q[0];
rz(2.334306) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(-0.92179006) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67547922) q[0];
sx q[0];
rz(-2.4371494) q[0];
sx q[0];
rz(2.9324313) q[0];
rz(-pi) q[1];
rz(-2.3085262) q[2];
sx q[2];
rz(-1.3563915) q[2];
sx q[2];
rz(3.0467141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39107716) q[1];
sx q[1];
rz(-2.8120496) q[1];
sx q[1];
rz(-0.68117546) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7796302) q[3];
sx q[3];
rz(-2.4224671) q[3];
sx q[3];
rz(1.5428839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.90199295) q[2];
sx q[2];
rz(-2.5805876) q[2];
sx q[2];
rz(1.1716589) q[2];
rz(-1.7840067) q[3];
sx q[3];
rz(-1.6849018) q[3];
sx q[3];
rz(-1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25306025) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(-1.6660447) q[0];
rz(-1.5400003) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(2.4618861) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25526991) q[0];
sx q[0];
rz(-2.5810044) q[0];
sx q[0];
rz(0.38829304) q[0];
rz(1.3557415) q[2];
sx q[2];
rz(-0.70677033) q[2];
sx q[2];
rz(0.099345318) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6913773) q[1];
sx q[1];
rz(-0.9141578) q[1];
sx q[1];
rz(-2.0037829) q[1];
rz(-pi) q[2];
rz(2.1867832) q[3];
sx q[3];
rz(-0.98815742) q[3];
sx q[3];
rz(-0.59347502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1264964) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(-2.2311907) q[2];
rz(0.67534584) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(-0.85062406) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0891721) q[0];
sx q[0];
rz(-1.2344673) q[0];
sx q[0];
rz(-2.4269379) q[0];
rz(0.71406281) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(-1.9649327) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51417527) q[0];
sx q[0];
rz(-1.5222933) q[0];
sx q[0];
rz(-1.7343299) q[0];
rz(0.87877019) q[2];
sx q[2];
rz(-2.3792017) q[2];
sx q[2];
rz(2.9825485) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9709819) q[1];
sx q[1];
rz(-2.5787528) q[1];
sx q[1];
rz(-1.7112205) q[1];
rz(1.482974) q[3];
sx q[3];
rz(-2.0801968) q[3];
sx q[3];
rz(2.1109964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.24361336) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(-2.0142377) q[2];
rz(0.35774287) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7174299) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(-2.7453616) q[1];
sx q[1];
rz(-0.025066499) q[1];
sx q[1];
rz(0.33096663) q[1];
rz(1.8489807) q[2];
sx q[2];
rz(-0.85360151) q[2];
sx q[2];
rz(-3.1384946) q[2];
rz(0.022396537) q[3];
sx q[3];
rz(-0.35084421) q[3];
sx q[3];
rz(2.4435333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
