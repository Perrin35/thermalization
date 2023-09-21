OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(-1.3077868) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(4.0772822) q[1];
sx q[1];
rz(4.7128591) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4966272) q[0];
sx q[0];
rz(-2.6500406) q[0];
sx q[0];
rz(2.3636723) q[0];
x q[1];
rz(2.4351032) q[2];
sx q[2];
rz(-2.2405365) q[2];
sx q[2];
rz(-1.1342088) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.03304122) q[1];
sx q[1];
rz(-1.812495) q[1];
sx q[1];
rz(-2.7944195) q[1];
rz(-pi) q[2];
rz(2.0516112) q[3];
sx q[3];
rz(-0.40502031) q[3];
sx q[3];
rz(2.4570176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.87542614) q[2];
sx q[2];
rz(-0.29310075) q[2];
sx q[2];
rz(2.0092633) q[2];
rz(1.4663565) q[3];
sx q[3];
rz(-1.3365859) q[3];
sx q[3];
rz(-1.0124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19673008) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(0.18584132) q[0];
rz(-2.5813685) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(-2.9247608) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6235979) q[0];
sx q[0];
rz(-1.2835842) q[0];
sx q[0];
rz(2.2395796) q[0];
rz(0.88044135) q[2];
sx q[2];
rz(-0.67697064) q[2];
sx q[2];
rz(-1.134269) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9456957) q[1];
sx q[1];
rz(-0.90598124) q[1];
sx q[1];
rz(2.871454) q[1];
x q[2];
rz(-0.87644491) q[3];
sx q[3];
rz(-1.0507686) q[3];
sx q[3];
rz(-2.2976573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8314787) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(1.2878093) q[2];
rz(2.3790322) q[3];
sx q[3];
rz(-1.972714) q[3];
sx q[3];
rz(0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6771616) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(0.60423869) q[0];
rz(-1.3263946) q[1];
sx q[1];
rz(-1.7809968) q[1];
sx q[1];
rz(-2.2089829) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9661449) q[0];
sx q[0];
rz(-1.6245337) q[0];
sx q[0];
rz(1.2530112) q[0];
rz(-pi) q[1];
rz(2.9595397) q[2];
sx q[2];
rz(-1.3472054) q[2];
sx q[2];
rz(-1.455866) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31308094) q[1];
sx q[1];
rz(-0.57144895) q[1];
sx q[1];
rz(2.9213195) q[1];
x q[2];
rz(1.014939) q[3];
sx q[3];
rz(-1.9564637) q[3];
sx q[3];
rz(-0.25394299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9937667) q[2];
sx q[2];
rz(-1.0819165) q[2];
sx q[2];
rz(2.0489342) q[2];
rz(2.5993733) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7820691) q[0];
sx q[0];
rz(-3.0451267) q[0];
sx q[0];
rz(-2.6413667) q[0];
rz(0.80530986) q[1];
sx q[1];
rz(-1.1601245) q[1];
sx q[1];
rz(1.6436228) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7786176) q[0];
sx q[0];
rz(-2.551429) q[0];
sx q[0];
rz(2.1182563) q[0];
rz(0.056604071) q[2];
sx q[2];
rz(-1.3805693) q[2];
sx q[2];
rz(-0.20400001) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8343463) q[1];
sx q[1];
rz(-1.3386968) q[1];
sx q[1];
rz(0.26718617) q[1];
rz(0.46521503) q[3];
sx q[3];
rz(-1.9808931) q[3];
sx q[3];
rz(-2.0801534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3952289) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(0.70181075) q[2];
rz(0.83135215) q[3];
sx q[3];
rz(-0.96389198) q[3];
sx q[3];
rz(0.62197661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2410626) q[0];
sx q[0];
rz(-2.5456972) q[0];
sx q[0];
rz(2.3262614) q[0];
rz(1.5218081) q[1];
sx q[1];
rz(-0.83414572) q[1];
sx q[1];
rz(-2.0933847) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74131504) q[0];
sx q[0];
rz(-1.321723) q[0];
sx q[0];
rz(-1.2988017) q[0];
rz(-pi) q[1];
rz(3.1254966) q[2];
sx q[2];
rz(-0.65158366) q[2];
sx q[2];
rz(-0.96166699) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.198846) q[1];
sx q[1];
rz(-1.8837351) q[1];
sx q[1];
rz(1.320977) q[1];
rz(-pi) q[2];
rz(-0.95543315) q[3];
sx q[3];
rz(-1.9121998) q[3];
sx q[3];
rz(1.1036901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.52577019) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(1.099951) q[2];
rz(-0.82529092) q[3];
sx q[3];
rz(-2.0992978) q[3];
sx q[3];
rz(0.88551372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0681756) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(0.90240479) q[0];
rz(-1.0166608) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(-3.0117603) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79486217) q[0];
sx q[0];
rz(-2.4299893) q[0];
sx q[0];
rz(0.58563389) q[0];
rz(2.1075222) q[2];
sx q[2];
rz(-0.64646361) q[2];
sx q[2];
rz(1.5922286) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2643913) q[1];
sx q[1];
rz(-1.0069205) q[1];
sx q[1];
rz(-1.1370204) q[1];
rz(-pi) q[2];
x q[2];
rz(1.767166) q[3];
sx q[3];
rz(-1.0305627) q[3];
sx q[3];
rz(1.2353209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3123902) q[2];
sx q[2];
rz(-2.1924993) q[2];
sx q[2];
rz(2.9373346) q[2];
rz(-1.9355109) q[3];
sx q[3];
rz(-1.5217425) q[3];
sx q[3];
rz(2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181353) q[0];
sx q[0];
rz(-1.3292987) q[0];
sx q[0];
rz(1.4468505) q[0];
rz(1.2591259) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(0.68626219) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0670358) q[0];
sx q[0];
rz(-2.0445163) q[0];
sx q[0];
rz(-1.0780328) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83001901) q[2];
sx q[2];
rz(-2.0247211) q[2];
sx q[2];
rz(1.1656851) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4916617) q[1];
sx q[1];
rz(-1.6451391) q[1];
sx q[1];
rz(-1.7564303) q[1];
rz(-pi) q[2];
rz(-2.7191914) q[3];
sx q[3];
rz(-1.8654612) q[3];
sx q[3];
rz(-0.15448031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4454322) q[2];
sx q[2];
rz(-1.3779209) q[2];
sx q[2];
rz(3.1398204) q[2];
rz(0.56162515) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(1.6368438) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6034265) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(1.4461393) q[0];
rz(0.7810477) q[1];
sx q[1];
rz(-1.8361517) q[1];
sx q[1];
rz(1.6400281) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43270375) q[0];
sx q[0];
rz(-0.49312691) q[0];
sx q[0];
rz(-1.5296442) q[0];
x q[1];
rz(1.1649706) q[2];
sx q[2];
rz(-3.0743982) q[2];
sx q[2];
rz(-0.42025987) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1135243) q[1];
sx q[1];
rz(-1.6599732) q[1];
sx q[1];
rz(2.813617) q[1];
rz(-pi) q[2];
rz(2.4092259) q[3];
sx q[3];
rz(-0.77419705) q[3];
sx q[3];
rz(-0.68294169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7897196) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(1.3191351) q[2];
rz(-1.2119279) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33655745) q[0];
sx q[0];
rz(-2.5890077) q[0];
sx q[0];
rz(-1.9375027) q[0];
rz(0.38326344) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(2.7899172) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4597804) q[0];
sx q[0];
rz(-2.5373055) q[0];
sx q[0];
rz(2.2629645) q[0];
x q[1];
rz(1.8750538) q[2];
sx q[2];
rz(-0.90196246) q[2];
sx q[2];
rz(-2.5077016) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5896776) q[1];
sx q[1];
rz(-0.81112408) q[1];
sx q[1];
rz(3.0685436) q[1];
rz(-pi) q[2];
rz(-2.6978108) q[3];
sx q[3];
rz(-0.1212596) q[3];
sx q[3];
rz(1.4951984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.3433156) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(1.2822255) q[2];
rz(-1.4964237) q[3];
sx q[3];
rz(-1.6069501) q[3];
sx q[3];
rz(-1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6431817) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(0.19432755) q[0];
rz(1.0378029) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(-1.0338354) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4612504) q[0];
sx q[0];
rz(-2.470812) q[0];
sx q[0];
rz(-1.781342) q[0];
rz(-pi) q[1];
rz(-0.50844426) q[2];
sx q[2];
rz(-2.2061081) q[2];
sx q[2];
rz(-0.19257643) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0018113) q[1];
sx q[1];
rz(-1.0011295) q[1];
sx q[1];
rz(-3.0210178) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8909573) q[3];
sx q[3];
rz(-0.72892979) q[3];
sx q[3];
rz(2.7540516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0795435) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(-0.6357843) q[2];
rz(2.87129) q[3];
sx q[3];
rz(-2.342194) q[3];
sx q[3];
rz(1.6132145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4476267) q[0];
sx q[0];
rz(-1.3128558) q[0];
sx q[0];
rz(-2.0679612) q[0];
rz(-1.7059965) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(-1.5244665) q[2];
sx q[2];
rz(-2.5382858) q[2];
sx q[2];
rz(-0.45679191) q[2];
rz(0.070449645) q[3];
sx q[3];
rz(-1.9000713) q[3];
sx q[3];
rz(-2.6303359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
