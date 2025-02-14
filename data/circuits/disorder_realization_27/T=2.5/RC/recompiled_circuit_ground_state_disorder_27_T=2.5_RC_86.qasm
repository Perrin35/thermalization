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
rz(-0.27443767) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4532758) q[0];
sx q[0];
rz(-0.098628086) q[0];
sx q[0];
rz(-2.3345956) q[0];
rz(-1.1054613) q[2];
sx q[2];
rz(-1.0792152) q[2];
sx q[2];
rz(2.9228989) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5861862) q[1];
sx q[1];
rz(-1.6073391) q[1];
sx q[1];
rz(-0.010374109) q[1];
rz(-pi) q[2];
rz(1.7414344) q[3];
sx q[3];
rz(-1.4981761) q[3];
sx q[3];
rz(-2.448248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9258257) q[2];
sx q[2];
rz(-3.1302858) q[2];
sx q[2];
rz(-2.0634148) q[2];
rz(-0.82873851) q[3];
sx q[3];
rz(-1.6308035) q[3];
sx q[3];
rz(-0.80882788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0593798) q[0];
sx q[0];
rz(-1.2780715) q[0];
sx q[0];
rz(-1.7510121) q[0];
rz(1.4311721) q[1];
sx q[1];
rz(-0.0043967604) q[1];
sx q[1];
rz(1.4336525) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2947049) q[0];
sx q[0];
rz(-1.6754912) q[0];
sx q[0];
rz(-0.93331915) q[0];
rz(-pi) q[1];
rz(-1.6054691) q[2];
sx q[2];
rz(-2.6975687) q[2];
sx q[2];
rz(-1.498986) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8987309) q[1];
sx q[1];
rz(-0.0090829385) q[1];
sx q[1];
rz(1.6159978) q[1];
rz(-0.065304718) q[3];
sx q[3];
rz(-0.78867824) q[3];
sx q[3];
rz(0.95897934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37890515) q[2];
sx q[2];
rz(-1.5964369) q[2];
sx q[2];
rz(-1.5691266) q[2];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.222027) q[0];
sx q[0];
rz(-2.5313105) q[0];
sx q[0];
rz(0.56184226) q[0];
rz(-1.5737083) q[1];
sx q[1];
rz(-1.5627197) q[1];
sx q[1];
rz(0.017339658) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0837096) q[0];
sx q[0];
rz(-1.0824993) q[0];
sx q[0];
rz(2.6668307) q[0];
rz(-pi) q[1];
rz(0.88490136) q[2];
sx q[2];
rz(-1.9608616) q[2];
sx q[2];
rz(1.5221514) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9063532) q[1];
sx q[1];
rz(-1.5867375) q[1];
sx q[1];
rz(1.9332744) q[1];
rz(-pi) q[2];
x q[2];
rz(0.041009237) q[3];
sx q[3];
rz(-2.6309391) q[3];
sx q[3];
rz(-1.2488431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6277546) q[2];
sx q[2];
rz(-1.7226115) q[2];
sx q[2];
rz(0.55297744) q[2];
rz(-1.9364457) q[3];
sx q[3];
rz(-1.5670992) q[3];
sx q[3];
rz(-1.5945826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7503081) q[0];
sx q[0];
rz(-1.0306083) q[0];
sx q[0];
rz(1.3543825) q[0];
rz(1.3722108) q[1];
sx q[1];
rz(-3.1387699) q[1];
sx q[1];
rz(1.7599958) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041080495) q[0];
sx q[0];
rz(-1.1100475) q[0];
sx q[0];
rz(-2.6897088) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4468378) q[2];
sx q[2];
rz(-3.1379897) q[2];
sx q[2];
rz(2.3191593) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0160675) q[1];
sx q[1];
rz(-1.9125332) q[1];
sx q[1];
rz(2.2831232) q[1];
rz(-0.80251191) q[3];
sx q[3];
rz(-1.3617235) q[3];
sx q[3];
rz(-1.0788364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0733205) q[2];
sx q[2];
rz(-3.1219411) q[2];
sx q[2];
rz(1.9015296) q[2];
rz(-0.23010075) q[3];
sx q[3];
rz(-0.0041882526) q[3];
sx q[3];
rz(2.7219462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.087990046) q[0];
sx q[0];
rz(-0.78730655) q[0];
sx q[0];
rz(1.8863652) q[0];
rz(3.1322196) q[1];
sx q[1];
rz(-1.7724937) q[1];
sx q[1];
rz(-3.110041) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9242212) q[0];
sx q[0];
rz(-1.6243906) q[0];
sx q[0];
rz(-1.0384667) q[0];
rz(1.2818665) q[2];
sx q[2];
rz(-1.4865371) q[2];
sx q[2];
rz(0.97848985) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42571354) q[1];
sx q[1];
rz(-1.7919994) q[1];
sx q[1];
rz(-1.0278661) q[1];
x q[2];
rz(-0.95057861) q[3];
sx q[3];
rz(-2.0671386) q[3];
sx q[3];
rz(2.2994808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.3343398) q[2];
sx q[2];
rz(-3.1355317) q[2];
sx q[2];
rz(-2.3414211) q[2];
rz(-2.3470894) q[3];
sx q[3];
rz(-0.032363351) q[3];
sx q[3];
rz(-2.2728424) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.04258) q[0];
sx q[0];
rz(-0.96262162) q[0];
sx q[0];
rz(0.30607589) q[0];
rz(-1.4535657) q[2];
sx q[2];
rz(-2.3701326) q[2];
sx q[2];
rz(-1.3958803) q[2];
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
rz(-pi) q[2];
x q[2];
rz(-2.6088151) q[3];
sx q[3];
rz(-2.2190337) q[3];
sx q[3];
rz(-1.4545573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.738203) q[2];
sx q[2];
rz(-0.047813606) q[2];
sx q[2];
rz(1.8955463) q[2];
rz(-1.3561148) q[3];
sx q[3];
rz(-3.1062283) q[3];
sx q[3];
rz(-1.416392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4288915) q[0];
sx q[0];
rz(-2.3172947) q[0];
sx q[0];
rz(-1.4225381) q[0];
rz(-1.2369583) q[1];
sx q[1];
rz(-3.1045541) q[1];
sx q[1];
rz(-2.9388156) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1272498) q[0];
sx q[0];
rz(-0.83006751) q[0];
sx q[0];
rz(-1.1237445) q[0];
rz(-2.3788664) q[2];
sx q[2];
rz(-0.87050754) q[2];
sx q[2];
rz(1.8000079) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8824892) q[1];
sx q[1];
rz(-1.9111425) q[1];
sx q[1];
rz(1.629414) q[1];
x q[2];
rz(2.9843627) q[3];
sx q[3];
rz(-1.0649324) q[3];
sx q[3];
rz(-1.0613522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5286336) q[2];
sx q[2];
rz(-3.0411159) q[2];
sx q[2];
rz(-2.4670777) q[2];
rz(-1.7965192) q[3];
sx q[3];
rz(-0.14480545) q[3];
sx q[3];
rz(1.8176746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26789185) q[0];
sx q[0];
rz(-0.75650263) q[0];
sx q[0];
rz(-2.2896413) q[0];
rz(2.9497228) q[1];
sx q[1];
rz(-0.012902915) q[1];
sx q[1];
rz(0.26564863) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0276447) q[0];
sx q[0];
rz(-2.6922604) q[0];
sx q[0];
rz(0.41037257) q[0];
x q[1];
rz(-0.4531817) q[2];
sx q[2];
rz(-2.4754449) q[2];
sx q[2];
rz(2.0311126) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.78662465) q[1];
sx q[1];
rz(-1.6296128) q[1];
sx q[1];
rz(-3.0579964) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1962131) q[3];
sx q[3];
rz(-2.0959955) q[3];
sx q[3];
rz(-1.7708667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.77136451) q[2];
sx q[2];
rz(-3.0493272) q[2];
sx q[2];
rz(2.6293758) q[2];
rz(2.9825315) q[3];
sx q[3];
rz(-0.03511196) q[3];
sx q[3];
rz(-1.4353282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-2.2559741) q[0];
sx q[0];
rz(-1.6640478) q[0];
sx q[0];
rz(-1.0783827) q[0];
rz(1.501561) q[1];
sx q[1];
rz(-2.9402132) q[1];
sx q[1];
rz(-1.5578425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5659065) q[0];
sx q[0];
rz(-2.4543336) q[0];
sx q[0];
rz(0.024373011) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89084004) q[2];
sx q[2];
rz(-1.8906697) q[2];
sx q[2];
rz(-0.25242463) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.16405995) q[1];
sx q[1];
rz(-1.5741473) q[1];
sx q[1];
rz(-0.022025755) q[1];
rz(-1.1224062) q[3];
sx q[3];
rz(-2.4351154) q[3];
sx q[3];
rz(-1.1170596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.885159) q[2];
sx q[2];
rz(-3.1270471) q[2];
sx q[2];
rz(3.0719768) q[2];
rz(-0.066667892) q[3];
sx q[3];
rz(-1.0144517) q[3];
sx q[3];
rz(2.4623509) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2081864) q[0];
sx q[0];
rz(-1.8649768) q[0];
sx q[0];
rz(0.25190121) q[0];
rz(-1.6690286) q[1];
sx q[1];
rz(-2.9258969) q[1];
sx q[1];
rz(-3.0676945) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2011482) q[0];
sx q[0];
rz(-1.6415039) q[0];
sx q[0];
rz(-1.9364249) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5394707) q[2];
sx q[2];
rz(-1.5534473) q[2];
sx q[2];
rz(-0.67048873) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.628129) q[1];
sx q[1];
rz(-0.76490116) q[1];
sx q[1];
rz(2.7897863) q[1];
rz(0.20404317) q[3];
sx q[3];
rz(-1.4082485) q[3];
sx q[3];
rz(-2.7891087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.68022388) q[2];
sx q[2];
rz(-0.0071439925) q[2];
sx q[2];
rz(2.3738677) q[2];
rz(-1.7420306) q[3];
sx q[3];
rz(-0.00082409516) q[3];
sx q[3];
rz(-0.50423938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
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
rz(3.1172251) q[1];
sx q[1];
rz(-2.9822646) q[1];
sx q[1];
rz(0.23039625) q[1];
rz(-0.38300285) q[2];
sx q[2];
rz(-2.2097864) q[2];
sx q[2];
rz(-1.0351444) q[2];
rz(0.41718116) q[3];
sx q[3];
rz(-1.659703) q[3];
sx q[3];
rz(-0.27123527) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
