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
rz(1.525653) q[1];
sx q[1];
rz(4.6620044) q[1];
sx q[1];
rz(9.6992156) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0204791) q[0];
sx q[0];
rz(-1.502636) q[0];
sx q[0];
rz(-1.4994552) q[0];
x q[1];
rz(2.0361314) q[2];
sx q[2];
rz(-1.0792152) q[2];
sx q[2];
rz(2.9228989) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1573616) q[1];
sx q[1];
rz(-1.5604291) q[1];
sx q[1];
rz(-1.6073411) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7414344) q[3];
sx q[3];
rz(-1.4981761) q[3];
sx q[3];
rz(-0.69334465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9258257) q[2];
sx q[2];
rz(-0.01130686) q[2];
sx q[2];
rz(-1.0781778) q[2];
rz(-2.3128541) q[3];
sx q[3];
rz(-1.5107892) q[3];
sx q[3];
rz(-0.80882788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082212903) q[0];
sx q[0];
rz(-1.2780715) q[0];
sx q[0];
rz(1.3905806) q[0];
rz(-1.7104205) q[1];
sx q[1];
rz(-0.0043967604) q[1];
sx q[1];
rz(-1.7079401) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7253256) q[0];
sx q[0];
rz(-2.4967589) q[0];
sx q[0];
rz(-1.3960442) q[0];
rz(1.1270056) q[2];
sx q[2];
rz(-1.5559042) q[2];
sx q[2];
rz(3.1010951) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8588577) q[1];
sx q[1];
rz(-1.5712067) q[1];
sx q[1];
rz(1.57987) q[1];
rz(-pi) q[2];
rz(-1.6363899) q[3];
sx q[3];
rz(-2.3573313) q[3];
sx q[3];
rz(-2.0900871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7626875) q[2];
sx q[2];
rz(-1.5964369) q[2];
sx q[2];
rz(1.5691266) q[2];
rz(2.3804741) q[3];
sx q[3];
rz(-3.0877536) q[3];
sx q[3];
rz(1.9332473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.222027) q[0];
sx q[0];
rz(-2.5313105) q[0];
sx q[0];
rz(0.56184226) q[0];
rz(1.5737083) q[1];
sx q[1];
rz(-1.5627197) q[1];
sx q[1];
rz(3.124253) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27631381) q[0];
sx q[0];
rz(-1.9863578) q[0];
sx q[0];
rz(2.1091976) q[0];
rz(-pi) q[1];
rz(0.48835619) q[2];
sx q[2];
rz(-2.1966509) q[2];
sx q[2];
rz(2.7912315) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9063532) q[1];
sx q[1];
rz(-1.5548551) q[1];
sx q[1];
rz(1.9332744) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.041009237) q[3];
sx q[3];
rz(-0.51065356) q[3];
sx q[3];
rz(1.8927495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5138381) q[2];
sx q[2];
rz(-1.7226115) q[2];
sx q[2];
rz(-2.5886152) q[2];
rz(-1.9364457) q[3];
sx q[3];
rz(-1.5744934) q[3];
sx q[3];
rz(-1.5470101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39128458) q[0];
sx q[0];
rz(-2.1109844) q[0];
sx q[0];
rz(-1.7872101) q[0];
rz(1.3722108) q[1];
sx q[1];
rz(-3.1387699) q[1];
sx q[1];
rz(-1.3815968) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3171661) q[0];
sx q[0];
rz(-1.9726511) q[0];
sx q[0];
rz(2.0749932) q[0];
x q[1];
rz(-0.0027678522) q[2];
sx q[2];
rz(-1.5684897) q[2];
sx q[2];
rz(-0.053611343) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3051532) q[1];
sx q[1];
rz(-2.2341995) q[1];
sx q[1];
rz(0.43933489) q[1];
rz(-pi) q[2];
rz(0.80251191) q[3];
sx q[3];
rz(-1.3617235) q[3];
sx q[3];
rz(-2.0627562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0733205) q[2];
sx q[2];
rz(-0.019651532) q[2];
sx q[2];
rz(-1.9015296) q[2];
rz(0.23010075) q[3];
sx q[3];
rz(-0.0041882526) q[3];
sx q[3];
rz(-2.7219462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0536026) q[0];
sx q[0];
rz(-0.78730655) q[0];
sx q[0];
rz(-1.8863652) q[0];
rz(3.1322196) q[1];
sx q[1];
rz(-1.3690989) q[1];
sx q[1];
rz(-0.031551687) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.384969) q[0];
sx q[0];
rz(-2.1022804) q[0];
sx q[0];
rz(3.079412) q[0];
rz(-pi) q[1];
rz(1.8597262) q[2];
sx q[2];
rz(-1.6550555) q[2];
sx q[2];
rz(0.97848985) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1281368) q[1];
sx q[1];
rz(-2.0990879) q[1];
sx q[1];
rz(0.25685132) q[1];
rz(-2.3208951) q[3];
sx q[3];
rz(-2.3682171) q[3];
sx q[3];
rz(-0.14062961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8072529) q[2];
sx q[2];
rz(-0.0060609239) q[2];
sx q[2];
rz(-0.80017153) q[2];
rz(-0.79450327) q[3];
sx q[3];
rz(-3.1092293) q[3];
sx q[3];
rz(0.86875027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(-1.054258) q[0];
sx q[0];
rz(-2.9825409) q[0];
sx q[0];
rz(-1.6416838) q[0];
rz(0.17290393) q[1];
sx q[1];
rz(-3.0992295) q[1];
sx q[1];
rz(-0.049887966) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2931516) q[0];
sx q[0];
rz(-1.3209136) q[0];
sx q[0];
rz(2.2014653) q[0];
x q[1];
rz(2.3388145) q[2];
sx q[2];
rz(-1.652431) q[2];
sx q[2];
rz(0.25914524) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.0058635423) q[1];
sx q[1];
rz(-1.3220437) q[1];
sx q[1];
rz(2.2963524) q[1];
rz(-pi) q[2];
rz(2.2920554) q[3];
sx q[3];
rz(-1.153933) q[3];
sx q[3];
rz(0.45826926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40338966) q[2];
sx q[2];
rz(-3.093779) q[2];
sx q[2];
rz(1.8955463) q[2];
rz(1.7854779) q[3];
sx q[3];
rz(-0.035364371) q[3];
sx q[3];
rz(1.416392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7127011) q[0];
sx q[0];
rz(-2.3172947) q[0];
sx q[0];
rz(-1.4225381) q[0];
rz(-1.9046344) q[1];
sx q[1];
rz(-0.037038602) q[1];
sx q[1];
rz(0.2027771) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13064676) q[0];
sx q[0];
rz(-1.8955064) q[0];
sx q[0];
rz(-2.3492011) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76272623) q[2];
sx q[2];
rz(-0.87050754) q[2];
sx q[2];
rz(-1.8000079) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8494871) q[1];
sx q[1];
rz(-1.6260481) q[1];
sx q[1];
rz(2.8007052) q[1];
x q[2];
rz(1.8462798) q[3];
sx q[3];
rz(-0.52770319) q[3];
sx q[3];
rz(2.3964411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5286336) q[2];
sx q[2];
rz(-0.10047675) q[2];
sx q[2];
rz(-2.4670777) q[2];
rz(-1.7965192) q[3];
sx q[3];
rz(-0.14480545) q[3];
sx q[3];
rz(-1.323918) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8737008) q[0];
sx q[0];
rz(-0.75650263) q[0];
sx q[0];
rz(-0.85195136) q[0];
rz(-2.9497228) q[1];
sx q[1];
rz(-0.012902915) q[1];
sx q[1];
rz(2.875944) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11394792) q[0];
sx q[0];
rz(-0.4493323) q[0];
sx q[0];
rz(2.7312201) q[0];
x q[1];
rz(2.688411) q[2];
sx q[2];
rz(-0.66614775) q[2];
sx q[2];
rz(1.11048) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77924624) q[1];
sx q[1];
rz(-1.4873449) q[1];
sx q[1];
rz(1.6298184) q[1];
rz(-pi) q[2];
rz(0.56318483) q[3];
sx q[3];
rz(-2.5068589) q[3];
sx q[3];
rz(0.70574443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3702281) q[2];
sx q[2];
rz(-3.0493272) q[2];
sx q[2];
rz(-0.51221687) q[2];
rz(2.9825315) q[3];
sx q[3];
rz(-3.1064807) q[3];
sx q[3];
rz(1.4353282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2559741) q[0];
sx q[0];
rz(-1.6640478) q[0];
sx q[0];
rz(-2.0632099) q[0];
rz(-1.6400317) q[1];
sx q[1];
rz(-0.20137943) q[1];
sx q[1];
rz(1.5578425) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013951741) q[0];
sx q[0];
rz(-1.5553345) q[0];
sx q[0];
rz(2.4544793) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40270771) q[2];
sx q[2];
rz(-0.93120775) q[2];
sx q[2];
rz(-2.072203) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4068102) q[1];
sx q[1];
rz(-1.5487707) q[1];
sx q[1];
rz(1.5741482) q[1];
x q[2];
rz(0.91520941) q[3];
sx q[3];
rz(-1.2855227) q[3];
sx q[3];
rz(-2.3370621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.885159) q[2];
sx q[2];
rz(-3.1270471) q[2];
sx q[2];
rz(-3.0719768) q[2];
rz(3.0749248) q[3];
sx q[3];
rz(-1.0144517) q[3];
sx q[3];
rz(2.4623509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2081864) q[0];
sx q[0];
rz(-1.8649768) q[0];
sx q[0];
rz(-0.25190121) q[0];
rz(-1.4725641) q[1];
sx q[1];
rz(-0.21569574) q[1];
sx q[1];
rz(-3.0676945) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44786772) q[0];
sx q[0];
rz(-2.7694919) q[0];
sx q[0];
rz(-1.3752346) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5394707) q[2];
sx q[2];
rz(-1.5881453) q[2];
sx q[2];
rz(2.4711039) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.042797814) q[1];
sx q[1];
rz(-0.86314161) q[1];
sx q[1];
rz(1.8902197) q[1];
rz(0.20404317) q[3];
sx q[3];
rz(-1.7333441) q[3];
sx q[3];
rz(-0.35248392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.68022388) q[2];
sx q[2];
rz(-3.1344487) q[2];
sx q[2];
rz(-0.76772493) q[2];
rz(1.7420306) q[3];
sx q[3];
rz(-0.00082409516) q[3];
sx q[3];
rz(0.50423938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6298228) q[0];
sx q[0];
rz(-2.1586824) q[0];
sx q[0];
rz(-1.4221738) q[0];
rz(-0.024367532) q[1];
sx q[1];
rz(-2.9822646) q[1];
sx q[1];
rz(0.23039625) q[1];
rz(-2.7585898) q[2];
sx q[2];
rz(-0.93180626) q[2];
sx q[2];
rz(2.1064482) q[2];
rz(-0.41718116) q[3];
sx q[3];
rz(-1.4818896) q[3];
sx q[3];
rz(2.8703574) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
