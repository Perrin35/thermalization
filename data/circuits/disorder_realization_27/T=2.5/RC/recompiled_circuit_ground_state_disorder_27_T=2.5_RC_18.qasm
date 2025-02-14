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
rz(-0.023356181) q[0];
sx q[0];
rz(-2.2060702) q[0];
rz(-1.6159396) q[1];
sx q[1];
rz(-1.5204117) q[1];
sx q[1];
rz(2.867155) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0204791) q[0];
sx q[0];
rz(-1.6389567) q[0];
sx q[0];
rz(1.4994552) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69752946) q[2];
sx q[2];
rz(-2.4781057) q[2];
sx q[2];
rz(1.0349719) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5554065) q[1];
sx q[1];
rz(-1.6073391) q[1];
sx q[1];
rz(-0.010374109) q[1];
rz(-pi) q[2];
rz(-1.1660444) q[3];
sx q[3];
rz(-0.18530986) q[3];
sx q[3];
rz(2.6626793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9258257) q[2];
sx q[2];
rz(-3.1302858) q[2];
sx q[2];
rz(2.0634148) q[2];
rz(-0.82873851) q[3];
sx q[3];
rz(-1.6308035) q[3];
sx q[3];
rz(2.3327648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082212903) q[0];
sx q[0];
rz(-1.8635211) q[0];
sx q[0];
rz(1.7510121) q[0];
rz(-1.7104205) q[1];
sx q[1];
rz(-3.1371959) q[1];
sx q[1];
rz(-1.4336525) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7253256) q[0];
sx q[0];
rz(-0.64483374) q[0];
sx q[0];
rz(-1.7455484) q[0];
rz(-pi) q[1];
x q[1];
rz(0.016489224) q[2];
sx q[2];
rz(-1.1270583) q[2];
sx q[2];
rz(-1.6042142) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.28273496) q[1];
sx q[1];
rz(-1.5712067) q[1];
sx q[1];
rz(1.57987) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5052027) q[3];
sx q[3];
rz(-2.3573313) q[3];
sx q[3];
rz(-2.0900871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7626875) q[2];
sx q[2];
rz(-1.5964369) q[2];
sx q[2];
rz(1.572466) q[2];
rz(2.3804741) q[3];
sx q[3];
rz(-0.053839024) q[3];
sx q[3];
rz(-1.9332473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.222027) q[0];
sx q[0];
rz(-0.61028218) q[0];
sx q[0];
rz(-0.56184226) q[0];
rz(1.5678844) q[1];
sx q[1];
rz(-1.578873) q[1];
sx q[1];
rz(-0.017339658) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2527232) q[0];
sx q[0];
rz(-2.4742352) q[0];
sx q[0];
rz(0.86020893) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1465535) q[2];
sx q[2];
rz(-0.77313609) q[2];
sx q[2];
rz(-0.38627689) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8120809) q[1];
sx q[1];
rz(-1.9332262) q[1];
sx q[1];
rz(0.017048841) q[1];
x q[2];
rz(-1.5937599) q[3];
sx q[3];
rz(-2.0809789) q[3];
sx q[3];
rz(-1.2018454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5138381) q[2];
sx q[2];
rz(-1.7226115) q[2];
sx q[2];
rz(2.5886152) q[2];
rz(-1.9364457) q[3];
sx q[3];
rz(-1.5744934) q[3];
sx q[3];
rz(1.5945826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39128458) q[0];
sx q[0];
rz(-1.0306083) q[0];
sx q[0];
rz(-1.3543825) q[0];
rz(-1.7693819) q[1];
sx q[1];
rz(-0.0028227614) q[1];
sx q[1];
rz(-1.7599958) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3171661) q[0];
sx q[0];
rz(-1.1689416) q[0];
sx q[0];
rz(1.0665994) q[0];
rz(1.5684897) q[2];
sx q[2];
rz(-1.5680285) q[2];
sx q[2];
rz(1.6244013) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.18455566) q[1];
sx q[1];
rz(-0.77691764) q[1];
sx q[1];
rz(-2.0691815) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28691157) q[3];
sx q[3];
rz(-2.3182456) q[3];
sx q[3];
rz(-2.8475743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0733205) q[2];
sx q[2];
rz(-3.1219411) q[2];
sx q[2];
rz(1.2400631) q[2];
rz(-0.23010075) q[3];
sx q[3];
rz(-0.0041882526) q[3];
sx q[3];
rz(-0.41964644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(0.087990046) q[0];
sx q[0];
rz(-2.3542861) q[0];
sx q[0];
rz(1.8863652) q[0];
rz(-0.0093731006) q[1];
sx q[1];
rz(-1.7724937) q[1];
sx q[1];
rz(-3.110041) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2173715) q[0];
sx q[0];
rz(-1.6243906) q[0];
sx q[0];
rz(2.103126) q[0];
rz(-3.0537082) q[2];
sx q[2];
rz(-1.2829208) q[2];
sx q[2];
rz(-0.56729546) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1281368) q[1];
sx q[1];
rz(-1.0425048) q[1];
sx q[1];
rz(0.25685132) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.191014) q[3];
sx q[3];
rz(-2.0671386) q[3];
sx q[3];
rz(0.84211189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.3343398) q[2];
sx q[2];
rz(-3.1355317) q[2];
sx q[2];
rz(-0.80017153) q[2];
rz(2.3470894) q[3];
sx q[3];
rz(-0.032363351) q[3];
sx q[3];
rz(2.2728424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.054258) q[0];
sx q[0];
rz(-2.9825409) q[0];
sx q[0];
rz(1.4999088) q[0];
rz(0.17290393) q[1];
sx q[1];
rz(-0.042363107) q[1];
sx q[1];
rz(0.049887966) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5373851) q[0];
sx q[0];
rz(-2.4695463) q[0];
sx q[0];
rz(1.9792569) q[0];
rz(3.028333) q[2];
sx q[2];
rz(-2.335603) q[2];
sx q[2];
rz(-1.2330556) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3059002) q[1];
sx q[1];
rz(-0.75961194) q[1];
sx q[1];
rz(-1.2051969) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6088151) q[3];
sx q[3];
rz(-0.92255892) q[3];
sx q[3];
rz(1.6870354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40338966) q[2];
sx q[2];
rz(-0.047813606) q[2];
sx q[2];
rz(-1.2460463) q[2];
rz(-1.3561148) q[3];
sx q[3];
rz(-0.035364371) q[3];
sx q[3];
rz(1.416392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7127011) q[0];
sx q[0];
rz(-2.3172947) q[0];
sx q[0];
rz(-1.7190546) q[0];
rz(1.9046344) q[1];
sx q[1];
rz(-0.037038602) q[1];
sx q[1];
rz(2.9388156) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13064676) q[0];
sx q[0];
rz(-1.8955064) q[0];
sx q[0];
rz(0.79239158) q[0];
rz(-0.70901633) q[2];
sx q[2];
rz(-1.0142027) q[2];
sx q[2];
rz(-2.8188044) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29210556) q[1];
sx q[1];
rz(-1.6260481) q[1];
sx q[1];
rz(2.8007052) q[1];
rz(-pi) q[2];
rz(-1.2953128) q[3];
sx q[3];
rz(-0.52770319) q[3];
sx q[3];
rz(-0.74515158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5286336) q[2];
sx q[2];
rz(-0.10047675) q[2];
sx q[2];
rz(-0.67451492) q[2];
rz(1.3450735) q[3];
sx q[3];
rz(-2.9967872) q[3];
sx q[3];
rz(1.323918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26789185) q[0];
sx q[0];
rz(-2.38509) q[0];
sx q[0];
rz(-0.85195136) q[0];
rz(-2.9497228) q[1];
sx q[1];
rz(-0.012902915) q[1];
sx q[1];
rz(-0.26564863) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.083355) q[0];
sx q[0];
rz(-1.3966271) q[0];
sx q[0];
rz(2.7252498) q[0];
rz(-0.61518367) q[2];
sx q[2];
rz(-1.2968204) q[2];
sx q[2];
rz(-0.09466234) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.77924624) q[1];
sx q[1];
rz(-1.6542477) q[1];
sx q[1];
rz(1.6298184) q[1];
x q[2];
rz(2.5846768) q[3];
sx q[3];
rz(-1.2486826) q[3];
sx q[3];
rz(2.7469001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3702281) q[2];
sx q[2];
rz(-0.092265487) q[2];
sx q[2];
rz(2.6293758) q[2];
rz(2.9825315) q[3];
sx q[3];
rz(-0.03511196) q[3];
sx q[3];
rz(1.7062645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
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
rz(-2.9402132) q[1];
sx q[1];
rz(1.5578425) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5974332) q[0];
sx q[0];
rz(-0.88378105) q[0];
sx q[0];
rz(-1.5907955) q[0];
rz(-pi) q[1];
rz(2.2507526) q[2];
sx q[2];
rz(-1.250923) q[2];
sx q[2];
rz(2.889168) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4068102) q[1];
sx q[1];
rz(-1.592822) q[1];
sx q[1];
rz(1.5741482) q[1];
rz(-pi) q[2];
rz(-1.1224062) q[3];
sx q[3];
rz(-0.70647722) q[3];
sx q[3];
rz(1.1170596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.885159) q[2];
sx q[2];
rz(-3.1270471) q[2];
sx q[2];
rz(0.069615901) q[2];
rz(0.066667892) q[3];
sx q[3];
rz(-1.0144517) q[3];
sx q[3];
rz(0.67924172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2081864) q[0];
sx q[0];
rz(-1.2766159) q[0];
sx q[0];
rz(2.8896914) q[0];
rz(-1.4725641) q[1];
sx q[1];
rz(-2.9258969) q[1];
sx q[1];
rz(3.0676945) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94044444) q[0];
sx q[0];
rz(-1.6415039) q[0];
sx q[0];
rz(-1.9364249) q[0];
x q[1];
rz(2.0766856) q[2];
sx q[2];
rz(-3.1057851) q[2];
sx q[2];
rz(-2.7469025) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.042797814) q[1];
sx q[1];
rz(-0.86314161) q[1];
sx q[1];
rz(1.2513729) q[1];
x q[2];
rz(-2.4611887) q[3];
sx q[3];
rz(-0.26016737) q[3];
sx q[3];
rz(-0.55458595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.68022388) q[2];
sx q[2];
rz(-3.1344487) q[2];
sx q[2];
rz(2.3738677) q[2];
rz(-1.7420306) q[3];
sx q[3];
rz(-0.00082409516) q[3];
sx q[3];
rz(2.6373533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6298228) q[0];
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
rz(-1.6679933) q[3];
sx q[3];
rz(-1.9862277) q[3];
sx q[3];
rz(1.338892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
