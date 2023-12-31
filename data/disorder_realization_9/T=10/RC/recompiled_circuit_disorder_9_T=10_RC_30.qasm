OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7115241) q[0];
sx q[0];
rz(-0.067458955) q[0];
sx q[0];
rz(0.67396069) q[0];
rz(-0.31710467) q[1];
sx q[1];
rz(-1.6333406) q[1];
sx q[1];
rz(2.34692) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793517) q[0];
sx q[0];
rz(-1.9277713) q[0];
sx q[0];
rz(0.39874052) q[0];
x q[1];
rz(-1.4235714) q[2];
sx q[2];
rz(-1.0861673) q[2];
sx q[2];
rz(-1.1046315) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.049584576) q[1];
sx q[1];
rz(-1.2283748) q[1];
sx q[1];
rz(0.53978668) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0029293) q[3];
sx q[3];
rz(-2.0141467) q[3];
sx q[3];
rz(-1.5821622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.47544605) q[2];
sx q[2];
rz(-2.011516) q[2];
sx q[2];
rz(1.8165992) q[2];
rz(2.825286) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(-2.2345208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9544202) q[0];
sx q[0];
rz(-0.81048727) q[0];
sx q[0];
rz(2.1818838) q[0];
rz(1.487544) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(-0.62746343) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7146485) q[0];
sx q[0];
rz(-1.5639925) q[0];
sx q[0];
rz(-3.1338918) q[0];
x q[1];
rz(2.9299514) q[2];
sx q[2];
rz(-1.5923119) q[2];
sx q[2];
rz(-3.0880398) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5261425) q[1];
sx q[1];
rz(-0.28029385) q[1];
sx q[1];
rz(2.1089094) q[1];
rz(1.596132) q[3];
sx q[3];
rz(-1.5583894) q[3];
sx q[3];
rz(-0.62904639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54005694) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(2.1616948) q[2];
rz(1.4510138) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(-1.6987945) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56101218) q[0];
sx q[0];
rz(-3.09364) q[0];
sx q[0];
rz(1.3797492) q[0];
rz(2.9648932) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(-1.7938991) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028570024) q[0];
sx q[0];
rz(-1.1060113) q[0];
sx q[0];
rz(-1.8405429) q[0];
x q[1];
rz(-1.5930487) q[2];
sx q[2];
rz(-1.77102) q[2];
sx q[2];
rz(-0.74860886) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29979953) q[1];
sx q[1];
rz(-1.449297) q[1];
sx q[1];
rz(-1.4701162) q[1];
rz(-1.1191145) q[3];
sx q[3];
rz(-2.936755) q[3];
sx q[3];
rz(2.0126437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0195062) q[2];
sx q[2];
rz(-0.85947376) q[2];
sx q[2];
rz(1.1536095) q[2];
rz(0.38976088) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(0.78891689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6692114) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(-0.56458449) q[0];
rz(-3.0584884) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(-0.21534236) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9504844) q[0];
sx q[0];
rz(-1.2380301) q[0];
sx q[0];
rz(-0.82975181) q[0];
x q[1];
rz(1.8183069) q[2];
sx q[2];
rz(-2.5396721) q[2];
sx q[2];
rz(1.9844696) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.54289651) q[1];
sx q[1];
rz(-0.22194949) q[1];
sx q[1];
rz(2.0680244) q[1];
x q[2];
rz(-0.54939778) q[3];
sx q[3];
rz(-1.1664207) q[3];
sx q[3];
rz(-2.3466563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5061491) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(-2.5396458) q[2];
rz(2.4510032) q[3];
sx q[3];
rz(-2.3875321) q[3];
sx q[3];
rz(-2.5155892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(0.85201207) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(-0.96486282) q[0];
rz(3.124974) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(2.4749277) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4283838) q[0];
sx q[0];
rz(-0.41233006) q[0];
sx q[0];
rz(1.8250188) q[0];
x q[1];
rz(-0.35258106) q[2];
sx q[2];
rz(-1.0513628) q[2];
sx q[2];
rz(0.40371603) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.40068377) q[1];
sx q[1];
rz(-2.2212257) q[1];
sx q[1];
rz(1.9294192) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2090962) q[3];
sx q[3];
rz(-0.36358788) q[3];
sx q[3];
rz(2.5225085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.91606402) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(1.8886245) q[2];
rz(-1.3145674) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(1.130828) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0261633) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(1.3622267) q[0];
rz(1.399614) q[1];
sx q[1];
rz(-1.5796966) q[1];
sx q[1];
rz(-1.6113575) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8996457) q[0];
sx q[0];
rz(-0.42662222) q[0];
sx q[0];
rz(-2.9645779) q[0];
rz(-pi) q[1];
rz(0.43892626) q[2];
sx q[2];
rz(-1.8669087) q[2];
sx q[2];
rz(-0.93026464) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.02009) q[1];
sx q[1];
rz(-1.1768747) q[1];
sx q[1];
rz(-2.9083088) q[1];
rz(3.0909782) q[3];
sx q[3];
rz(-1.317987) q[3];
sx q[3];
rz(-1.8265754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7851012) q[2];
sx q[2];
rz(-1.0078398) q[2];
sx q[2];
rz(3.1121758) q[2];
rz(-1.4908837) q[3];
sx q[3];
rz(-2.1865032) q[3];
sx q[3];
rz(1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4239663) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(0.18280612) q[0];
rz(2.706066) q[1];
sx q[1];
rz(-2.2874449) q[1];
sx q[1];
rz(-2.6307154) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85833997) q[0];
sx q[0];
rz(-2.8562299) q[0];
sx q[0];
rz(2.7417389) q[0];
x q[1];
rz(-2.3235544) q[2];
sx q[2];
rz(-0.99070264) q[2];
sx q[2];
rz(1.8028508) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.43607298) q[1];
sx q[1];
rz(-1.9680068) q[1];
sx q[1];
rz(-0.7049837) q[1];
rz(-pi) q[2];
rz(-2.0914518) q[3];
sx q[3];
rz(-0.58590349) q[3];
sx q[3];
rz(1.0019765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9110979) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(-2.5210209) q[2];
rz(2.7190322) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34185103) q[0];
sx q[0];
rz(-2.7809679) q[0];
sx q[0];
rz(2.0525232) q[0];
rz(1.0808806) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(0.63527766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1875293) q[0];
sx q[0];
rz(-0.39217338) q[0];
sx q[0];
rz(1.1568882) q[0];
rz(-pi) q[1];
rz(1.8842472) q[2];
sx q[2];
rz(-2.6230276) q[2];
sx q[2];
rz(2.5093362) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0465225) q[1];
sx q[1];
rz(-0.97192837) q[1];
sx q[1];
rz(-1.2477161) q[1];
rz(-1.2653207) q[3];
sx q[3];
rz(-2.5442122) q[3];
sx q[3];
rz(1.3437831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.13715956) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(3.1414202) q[2];
rz(-2.4553283) q[3];
sx q[3];
rz(-0.72824794) q[3];
sx q[3];
rz(0.85787684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2962608) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(2.4849179) q[0];
rz(-0.36390057) q[1];
sx q[1];
rz(-2.6952422) q[1];
sx q[1];
rz(2.231853) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1247524) q[0];
sx q[0];
rz(-0.42718712) q[0];
sx q[0];
rz(0.70144002) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34809525) q[2];
sx q[2];
rz(-1.7762134) q[2];
sx q[2];
rz(1.5945438) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.64987953) q[1];
sx q[1];
rz(-2.8976106) q[1];
sx q[1];
rz(-2.6073543) q[1];
rz(-pi) q[2];
rz(-1.2551366) q[3];
sx q[3];
rz(-1.3254032) q[3];
sx q[3];
rz(-2.0592225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3108814) q[2];
sx q[2];
rz(-0.41822663) q[2];
sx q[2];
rz(0.48661423) q[2];
rz(0.1085554) q[3];
sx q[3];
rz(-2.4311549) q[3];
sx q[3];
rz(1.2238812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20031032) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(0.51668984) q[0];
rz(1.5962881) q[1];
sx q[1];
rz(-2.089112) q[1];
sx q[1];
rz(0.89458481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0827732) q[0];
sx q[0];
rz(-1.6130157) q[0];
sx q[0];
rz(-2.204338) q[0];
rz(-pi) q[1];
rz(-1.1012494) q[2];
sx q[2];
rz(-1.38648) q[2];
sx q[2];
rz(-1.4354401) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91720944) q[1];
sx q[1];
rz(-1.6341126) q[1];
sx q[1];
rz(1.9320095) q[1];
rz(2.9856332) q[3];
sx q[3];
rz(-2.2150063) q[3];
sx q[3];
rz(-2.9498266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.14780012) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(-1.0220035) q[2];
rz(-1.3379124) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5932896) q[0];
sx q[0];
rz(-1.6175445) q[0];
sx q[0];
rz(2.2486726) q[0];
rz(-2.100636) q[1];
sx q[1];
rz(-0.092408471) q[1];
sx q[1];
rz(1.668781) q[1];
rz(-0.9085761) q[2];
sx q[2];
rz(-2.3618345) q[2];
sx q[2];
rz(0.97710412) q[2];
rz(-1.8716807) q[3];
sx q[3];
rz(-1.5242929) q[3];
sx q[3];
rz(2.268189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
