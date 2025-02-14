OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55750027) q[0];
sx q[0];
rz(-3.1182365) q[0];
sx q[0];
rz(-0.93552247) q[0];
rz(-1.6159396) q[1];
sx q[1];
rz(-1.5204117) q[1];
sx q[1];
rz(2.867155) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4532758) q[0];
sx q[0];
rz(-0.098628086) q[0];
sx q[0];
rz(-2.3345956) q[0];
x q[1];
rz(-0.53977592) q[2];
sx q[2];
rz(-1.1640884) q[2];
sx q[2];
rz(1.1193898) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9842311) q[1];
sx q[1];
rz(-1.5604291) q[1];
sx q[1];
rz(1.6073411) q[1];
rz(-pi) q[2];
rz(-1.4001582) q[3];
sx q[3];
rz(-1.4981761) q[3];
sx q[3];
rz(0.69334465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2157669) q[2];
sx q[2];
rz(-3.1302858) q[2];
sx q[2];
rz(1.0781778) q[2];
rz(-0.82873851) q[3];
sx q[3];
rz(-1.6308035) q[3];
sx q[3];
rz(2.3327648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082212903) q[0];
sx q[0];
rz(-1.2780715) q[0];
sx q[0];
rz(1.3905806) q[0];
rz(-1.4311721) q[1];
sx q[1];
rz(-0.0043967604) q[1];
sx q[1];
rz(1.7079401) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7253256) q[0];
sx q[0];
rz(-0.64483374) q[0];
sx q[0];
rz(1.7455484) q[0];
rz(-pi) q[1];
rz(-3.1251034) q[2];
sx q[2];
rz(-2.0145344) q[2];
sx q[2];
rz(1.6042142) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8588577) q[1];
sx q[1];
rz(-1.5712067) q[1];
sx q[1];
rz(-1.5617227) q[1];
rz(-3.0762879) q[3];
sx q[3];
rz(-2.3529144) q[3];
sx q[3];
rz(0.95897934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37890515) q[2];
sx q[2];
rz(-1.5451558) q[2];
sx q[2];
rz(1.5691266) q[2];
rz(2.3804741) q[3];
sx q[3];
rz(-0.053839024) q[3];
sx q[3];
rz(1.2083453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.222027) q[0];
sx q[0];
rz(-2.5313105) q[0];
sx q[0];
rz(0.56184226) q[0];
rz(1.5678844) q[1];
sx q[1];
rz(-1.5627197) q[1];
sx q[1];
rz(0.017339658) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0578831) q[0];
sx q[0];
rz(-2.0590933) q[0];
sx q[0];
rz(-2.6668307) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6532365) q[2];
sx q[2];
rz(-0.94494176) q[2];
sx q[2];
rz(-0.35036119) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3775656) q[1];
sx q[1];
rz(-2.7787797) q[1];
sx q[1];
rz(-1.6157263) q[1];
x q[2];
rz(1.5478327) q[3];
sx q[3];
rz(-2.0809789) q[3];
sx q[3];
rz(-1.2018454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6277546) q[2];
sx q[2];
rz(-1.7226115) q[2];
sx q[2];
rz(0.55297744) q[2];
rz(1.2051469) q[3];
sx q[3];
rz(-1.5670992) q[3];
sx q[3];
rz(-1.5945826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39128458) q[0];
sx q[0];
rz(-2.1109844) q[0];
sx q[0];
rz(-1.7872101) q[0];
rz(-1.7693819) q[1];
sx q[1];
rz(-3.1387699) q[1];
sx q[1];
rz(-1.3815968) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1005122) q[0];
sx q[0];
rz(-2.0315451) q[0];
sx q[0];
rz(-0.45188388) q[0];
rz(-pi) q[1];
rz(-1.5731029) q[2];
sx q[2];
rz(-1.5735642) q[2];
sx q[2];
rz(-1.6244013) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.957037) q[1];
sx q[1];
rz(-0.77691764) q[1];
sx q[1];
rz(1.0724111) q[1];
x q[2];
rz(-2.8546811) q[3];
sx q[3];
rz(-0.82334703) q[3];
sx q[3];
rz(2.8475743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0682721) q[2];
sx q[2];
rz(-0.019651532) q[2];
sx q[2];
rz(-1.2400631) q[2];
rz(-2.9114919) q[3];
sx q[3];
rz(-0.0041882526) q[3];
sx q[3];
rz(0.41964644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0536026) q[0];
sx q[0];
rz(-0.78730655) q[0];
sx q[0];
rz(-1.8863652) q[0];
rz(-0.0093731006) q[1];
sx q[1];
rz(-1.3690989) q[1];
sx q[1];
rz(3.110041) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.878859) q[0];
sx q[0];
rz(-2.6068305) q[0];
sx q[0];
rz(1.6761024) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8589694) q[2];
sx q[2];
rz(-2.8409578) q[2];
sx q[2];
rz(-0.86821454) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6477709) q[1];
sx q[1];
rz(-0.58205172) q[1];
sx q[1];
rz(-1.9813367) q[1];
x q[2];
rz(-2.3208951) q[3];
sx q[3];
rz(-0.77337556) q[3];
sx q[3];
rz(0.14062961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8072529) q[2];
sx q[2];
rz(-3.1355317) q[2];
sx q[2];
rz(-2.3414211) q[2];
rz(0.79450327) q[3];
sx q[3];
rz(-3.1092293) q[3];
sx q[3];
rz(2.2728424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.054258) q[0];
sx q[0];
rz(-0.15905173) q[0];
sx q[0];
rz(1.6416838) q[0];
rz(-2.9686887) q[1];
sx q[1];
rz(-0.042363107) q[1];
sx q[1];
rz(0.049887966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099012696) q[0];
sx q[0];
rz(-2.178971) q[0];
sx q[0];
rz(0.30607589) q[0];
rz(-pi) q[1];
rz(-0.80277819) q[2];
sx q[2];
rz(-1.652431) q[2];
sx q[2];
rz(-2.8824474) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.0058635423) q[1];
sx q[1];
rz(-1.8195489) q[1];
sx q[1];
rz(-2.2963524) q[1];
rz(-0.53277758) q[3];
sx q[3];
rz(-2.2190337) q[3];
sx q[3];
rz(1.4545573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40338966) q[2];
sx q[2];
rz(-3.093779) q[2];
sx q[2];
rz(-1.8955463) q[2];
rz(1.3561148) q[3];
sx q[3];
rz(-3.1062283) q[3];
sx q[3];
rz(1.416392) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4288915) q[0];
sx q[0];
rz(-0.8242979) q[0];
sx q[0];
rz(1.4225381) q[0];
rz(-1.2369583) q[1];
sx q[1];
rz(-0.037038602) q[1];
sx q[1];
rz(2.9388156) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13064676) q[0];
sx q[0];
rz(-1.8955064) q[0];
sx q[0];
rz(0.79239158) q[0];
rz(-pi) q[1];
rz(-2.4325763) q[2];
sx q[2];
rz(-1.0142027) q[2];
sx q[2];
rz(2.8188044) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8494871) q[1];
sx q[1];
rz(-1.5155445) q[1];
sx q[1];
rz(2.8007052) q[1];
rz(-pi) q[2];
rz(-1.0596541) q[3];
sx q[3];
rz(-1.4333925) q[3];
sx q[3];
rz(-2.5554772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5286336) q[2];
sx q[2];
rz(-3.0411159) q[2];
sx q[2];
rz(-0.67451492) q[2];
rz(1.3450735) q[3];
sx q[3];
rz(-2.9967872) q[3];
sx q[3];
rz(-1.8176746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26789185) q[0];
sx q[0];
rz(-0.75650263) q[0];
sx q[0];
rz(2.2896413) q[0];
rz(-2.9497228) q[1];
sx q[1];
rz(-0.012902915) q[1];
sx q[1];
rz(-0.26564863) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11394792) q[0];
sx q[0];
rz(-2.6922604) q[0];
sx q[0];
rz(0.41037257) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9022361) q[2];
sx q[2];
rz(-2.1598615) q[2];
sx q[2];
rz(1.4765431) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7456671) q[1];
sx q[1];
rz(-0.10217459) q[1];
sx q[1];
rz(0.61421444) q[1];
rz(-pi) q[2];
rz(1.9453796) q[3];
sx q[3];
rz(-1.0455971) q[3];
sx q[3];
rz(1.370726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.77136451) q[2];
sx q[2];
rz(-0.092265487) q[2];
sx q[2];
rz(-0.51221687) q[2];
rz(-0.15906119) q[3];
sx q[3];
rz(-3.1064807) q[3];
sx q[3];
rz(1.4353282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2559741) q[0];
sx q[0];
rz(-1.6640478) q[0];
sx q[0];
rz(2.0632099) q[0];
rz(-1.501561) q[1];
sx q[1];
rz(-0.20137943) q[1];
sx q[1];
rz(1.5837502) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5756861) q[0];
sx q[0];
rz(-2.4543336) q[0];
sx q[0];
rz(-3.1172196) q[0];
x q[1];
rz(2.0556765) q[2];
sx q[2];
rz(-0.74046248) q[2];
sx q[2];
rz(1.6892576) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.16405995) q[1];
sx q[1];
rz(-1.5741473) q[1];
sx q[1];
rz(-0.022025755) q[1];
x q[2];
rz(0.91520941) q[3];
sx q[3];
rz(-1.85607) q[3];
sx q[3];
rz(2.3370621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.25643361) q[2];
sx q[2];
rz(-3.1270471) q[2];
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
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2081864) q[0];
sx q[0];
rz(-1.2766159) q[0];
sx q[0];
rz(-2.8896914) q[0];
rz(1.4725641) q[1];
sx q[1];
rz(-2.9258969) q[1];
sx q[1];
rz(-3.0676945) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65739261) q[0];
sx q[0];
rz(-1.2061242) q[0];
sx q[0];
rz(3.0658989) q[0];
rz(-pi) q[1];
x q[1];
rz(0.017357512) q[2];
sx q[2];
rz(-1.6021172) q[2];
sx q[2];
rz(-0.90085122) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0987948) q[1];
sx q[1];
rz(-0.86314161) q[1];
sx q[1];
rz(-1.8902197) q[1];
rz(-pi) q[2];
rz(1.4048673) q[3];
sx q[3];
rz(-1.36948) q[3];
sx q[3];
rz(1.8898058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4613688) q[2];
sx q[2];
rz(-3.1344487) q[2];
sx q[2];
rz(2.3738677) q[2];
rz(1.399562) q[3];
sx q[3];
rz(-0.00082409516) q[3];
sx q[3];
rz(2.6373533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6298228) q[0];
sx q[0];
rz(-2.1586824) q[0];
sx q[0];
rz(-1.4221738) q[0];
rz(3.1172251) q[1];
sx q[1];
rz(-2.9822646) q[1];
sx q[1];
rz(0.23039625) q[1];
rz(-1.1047614) q[2];
sx q[2];
rz(-2.4106422) q[2];
sx q[2];
rz(1.5123083) q[2];
rz(1.4735994) q[3];
sx q[3];
rz(-1.9862277) q[3];
sx q[3];
rz(1.338892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
