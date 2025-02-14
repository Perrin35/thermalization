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
rz(-0.29210583) q[0];
sx q[0];
rz(-0.6120683) q[0];
sx q[0];
rz(-3.0906313) q[0];
rz(-1.851097) q[1];
sx q[1];
rz(-1.8654856) q[1];
sx q[1];
rz(-2.2955503) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.030662) q[0];
sx q[0];
rz(-1.4101068) q[0];
sx q[0];
rz(2.7881505) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8522367) q[2];
sx q[2];
rz(-1.932992) q[2];
sx q[2];
rz(1.8487428) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5436094) q[1];
sx q[1];
rz(-0.56386891) q[1];
sx q[1];
rz(-2.6257674) q[1];
rz(-pi) q[2];
rz(-1.5624763) q[3];
sx q[3];
rz(-1.6826311) q[3];
sx q[3];
rz(-2.73561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7840665) q[2];
sx q[2];
rz(-0.94782031) q[2];
sx q[2];
rz(0.86304647) q[2];
rz(2.0283902) q[3];
sx q[3];
rz(-2.2200255) q[3];
sx q[3];
rz(-2.34288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6031826) q[0];
sx q[0];
rz(-0.69284678) q[0];
sx q[0];
rz(0.26042724) q[0];
rz(2.8159091) q[1];
sx q[1];
rz(-0.98418701) q[1];
sx q[1];
rz(2.8338285) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4474898) q[0];
sx q[0];
rz(-1.8304951) q[0];
sx q[0];
rz(-2.3070241) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2529875) q[2];
sx q[2];
rz(-1.121064) q[2];
sx q[2];
rz(-2.1572402) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.766749) q[1];
sx q[1];
rz(-1.6628712) q[1];
sx q[1];
rz(1.8122775) q[1];
rz(-pi) q[2];
rz(2.2062088) q[3];
sx q[3];
rz(-1.5652802) q[3];
sx q[3];
rz(2.8817906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8327568) q[2];
sx q[2];
rz(-2.2726161) q[2];
sx q[2];
rz(-0.55574179) q[2];
rz(-2.6293788) q[3];
sx q[3];
rz(-0.54849505) q[3];
sx q[3];
rz(-0.51928025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25278768) q[0];
sx q[0];
rz(-0.35211173) q[0];
sx q[0];
rz(0.74932253) q[0];
rz(1.0961756) q[1];
sx q[1];
rz(-1.3744033) q[1];
sx q[1];
rz(0.38415092) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6644728) q[0];
sx q[0];
rz(-0.64711231) q[0];
sx q[0];
rz(2.6958334) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36142771) q[2];
sx q[2];
rz(-2.5088032) q[2];
sx q[2];
rz(0.91362) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2652581) q[1];
sx q[1];
rz(-1.6240082) q[1];
sx q[1];
rz(1.268178) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2189193) q[3];
sx q[3];
rz(-0.71563827) q[3];
sx q[3];
rz(-1.1688237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7595547) q[2];
sx q[2];
rz(-1.8363154) q[2];
sx q[2];
rz(1.6273512) q[2];
rz(0.73005992) q[3];
sx q[3];
rz(-0.74159139) q[3];
sx q[3];
rz(-2.7470284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.45383129) q[0];
sx q[0];
rz(-1.4241968) q[0];
sx q[0];
rz(-3.0815531) q[0];
rz(2.2278348) q[1];
sx q[1];
rz(-2.2303631) q[1];
sx q[1];
rz(-0.93889108) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70907043) q[0];
sx q[0];
rz(-0.071252206) q[0];
sx q[0];
rz(-0.48686102) q[0];
rz(-3.0591776) q[2];
sx q[2];
rz(-1.5396828) q[2];
sx q[2];
rz(-1.7331725) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42753427) q[1];
sx q[1];
rz(-2.6529544) q[1];
sx q[1];
rz(1.6227551) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4603953) q[3];
sx q[3];
rz(-1.0170611) q[3];
sx q[3];
rz(0.21986976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1303723) q[2];
sx q[2];
rz(-0.29272407) q[2];
sx q[2];
rz(2.4893153) q[2];
rz(1.8782714) q[3];
sx q[3];
rz(-1.5247366) q[3];
sx q[3];
rz(-1.1711082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071282722) q[0];
sx q[0];
rz(-1.2391397) q[0];
sx q[0];
rz(2.5791445) q[0];
rz(1.7930454) q[1];
sx q[1];
rz(-2.8548073) q[1];
sx q[1];
rz(-2.5602692) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5452257) q[0];
sx q[0];
rz(-2.2196182) q[0];
sx q[0];
rz(-2.5741386) q[0];
rz(-pi) q[1];
rz(-3.0688672) q[2];
sx q[2];
rz(-2.680034) q[2];
sx q[2];
rz(2.3978021) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6186468) q[1];
sx q[1];
rz(-0.60877555) q[1];
sx q[1];
rz(0.91880111) q[1];
rz(1.526889) q[3];
sx q[3];
rz(-0.9287408) q[3];
sx q[3];
rz(-1.8521483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.44744667) q[2];
sx q[2];
rz(-1.7116825) q[2];
sx q[2];
rz(-2.4549129) q[2];
rz(-0.40003362) q[3];
sx q[3];
rz(-1.8828078) q[3];
sx q[3];
rz(3.1299642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43468633) q[0];
sx q[0];
rz(-0.97395581) q[0];
sx q[0];
rz(2.6913225) q[0];
rz(-0.88092583) q[1];
sx q[1];
rz(-1.2342781) q[1];
sx q[1];
rz(-1.1164104) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6754782) q[0];
sx q[0];
rz(-0.7921392) q[0];
sx q[0];
rz(-1.303543) q[0];
x q[1];
rz(-1.3853778) q[2];
sx q[2];
rz(-1.2259019) q[2];
sx q[2];
rz(-1.138681) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19710572) q[1];
sx q[1];
rz(-2.53105) q[1];
sx q[1];
rz(-1.6050299) q[1];
rz(-pi) q[2];
rz(-2.0209337) q[3];
sx q[3];
rz(-1.880247) q[3];
sx q[3];
rz(1.3510973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.474596) q[2];
sx q[2];
rz(-2.907739) q[2];
sx q[2];
rz(-0.65143603) q[2];
rz(2.3598119) q[3];
sx q[3];
rz(-1.4835446) q[3];
sx q[3];
rz(2.206291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.5930138) q[0];
sx q[0];
rz(-0.23603708) q[0];
sx q[0];
rz(2.6759942) q[0];
rz(-1.9904526) q[1];
sx q[1];
rz(-1.8955889) q[1];
sx q[1];
rz(-1.8021072) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5975153) q[0];
sx q[0];
rz(-1.9751514) q[0];
sx q[0];
rz(-0.57283516) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8546245) q[2];
sx q[2];
rz(-3.0833088) q[2];
sx q[2];
rz(2.5789277) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0995746) q[1];
sx q[1];
rz(-2.9093643) q[1];
sx q[1];
rz(1.2378511) q[1];
x q[2];
rz(0.42134704) q[3];
sx q[3];
rz(-0.94910062) q[3];
sx q[3];
rz(-2.5558215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0307978) q[2];
sx q[2];
rz(-1.75533) q[2];
sx q[2];
rz(0.01290713) q[2];
rz(0.13478336) q[3];
sx q[3];
rz(-1.8811992) q[3];
sx q[3];
rz(-0.24136647) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4607234) q[0];
sx q[0];
rz(-3.0496106) q[0];
sx q[0];
rz(1.1750093) q[0];
rz(-2.3911632) q[1];
sx q[1];
rz(-1.78396) q[1];
sx q[1];
rz(-0.39355412) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65312306) q[0];
sx q[0];
rz(-0.61310378) q[0];
sx q[0];
rz(1.2638919) q[0];
rz(2.9904994) q[2];
sx q[2];
rz(-0.6616306) q[2];
sx q[2];
rz(-2.2349295) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.30545772) q[1];
sx q[1];
rz(-1.6806446) q[1];
sx q[1];
rz(-2.0019462) q[1];
rz(-pi) q[2];
rz(-0.92994382) q[3];
sx q[3];
rz(-0.612878) q[3];
sx q[3];
rz(1.5604492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63341081) q[2];
sx q[2];
rz(-1.1339302) q[2];
sx q[2];
rz(-2.9122162) q[2];
rz(3.0502099) q[3];
sx q[3];
rz(-1.4437851) q[3];
sx q[3];
rz(-0.55575931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2185739) q[0];
sx q[0];
rz(-0.78892437) q[0];
sx q[0];
rz(2.5564585) q[0];
rz(1.8642037) q[1];
sx q[1];
rz(-1.9265415) q[1];
sx q[1];
rz(2.9871984) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0994249) q[0];
sx q[0];
rz(-1.6800361) q[0];
sx q[0];
rz(2.7043998) q[0];
rz(-0.4265151) q[2];
sx q[2];
rz(-1.2517445) q[2];
sx q[2];
rz(1.691455) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.33340633) q[1];
sx q[1];
rz(-1.3097313) q[1];
sx q[1];
rz(-0.33269791) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.102907) q[3];
sx q[3];
rz(-1.7280735) q[3];
sx q[3];
rz(-3.0426971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.020891) q[2];
sx q[2];
rz(-0.88940826) q[2];
sx q[2];
rz(-0.54023877) q[2];
rz(0.39515105) q[3];
sx q[3];
rz(-0.65427798) q[3];
sx q[3];
rz(1.7925709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9304792) q[0];
sx q[0];
rz(-1.8607288) q[0];
sx q[0];
rz(1.5323918) q[0];
rz(-0.92597517) q[1];
sx q[1];
rz(-1.7151058) q[1];
sx q[1];
rz(1.6645974) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2944082) q[0];
sx q[0];
rz(-1.5700746) q[0];
sx q[0];
rz(-2.4154354) q[0];
rz(-pi) q[1];
rz(2.8468644) q[2];
sx q[2];
rz(-1.9325614) q[2];
sx q[2];
rz(2.8579762) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.6546911) q[1];
sx q[1];
rz(-2.4729252) q[1];
sx q[1];
rz(0.99744852) q[1];
rz(-0.98201237) q[3];
sx q[3];
rz(-1.9470306) q[3];
sx q[3];
rz(-2.328095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8605139) q[2];
sx q[2];
rz(-1.3382567) q[2];
sx q[2];
rz(0.64858428) q[2];
rz(2.4159238) q[3];
sx q[3];
rz(-2.9026493) q[3];
sx q[3];
rz(1.3919977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7657179) q[0];
sx q[0];
rz(-0.97828843) q[0];
sx q[0];
rz(0.82213415) q[0];
rz(1.0646461) q[1];
sx q[1];
rz(-1.6864265) q[1];
sx q[1];
rz(2.0815157) q[1];
rz(0.63659414) q[2];
sx q[2];
rz(-1.4298212) q[2];
sx q[2];
rz(0.22967568) q[2];
rz(-2.3271046) q[3];
sx q[3];
rz(-1.9662439) q[3];
sx q[3];
rz(1.6821485) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
