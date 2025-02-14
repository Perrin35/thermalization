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
rz(1.1887551) q[0];
sx q[0];
rz(-2.3030757) q[0];
sx q[0];
rz(1.4128348) q[0];
rz(0.38263327) q[1];
sx q[1];
rz(1.1525947) q[1];
sx q[1];
rz(11.424185) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6467616) q[0];
sx q[0];
rz(-1.9660304) q[0];
sx q[0];
rz(3.0656205) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0549802) q[2];
sx q[2];
rz(-1.118045) q[2];
sx q[2];
rz(-2.7198359) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.176291) q[1];
sx q[1];
rz(-1.3421196) q[1];
sx q[1];
rz(-2.2139174) q[1];
rz(-pi) q[2];
rz(2.258681) q[3];
sx q[3];
rz(-1.4174682) q[3];
sx q[3];
rz(0.49082498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6483868) q[2];
sx q[2];
rz(-1.1492665) q[2];
sx q[2];
rz(-0.0041848103) q[2];
rz(-2.3936213) q[3];
sx q[3];
rz(-2.3145521) q[3];
sx q[3];
rz(2.9048257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10906049) q[0];
sx q[0];
rz(-0.91962686) q[0];
sx q[0];
rz(-3.076886) q[0];
rz(1.0626571) q[1];
sx q[1];
rz(-0.92667842) q[1];
sx q[1];
rz(1.0437171) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4499265) q[0];
sx q[0];
rz(-2.0332893) q[0];
sx q[0];
rz(1.3426379) q[0];
rz(-3.114568) q[2];
sx q[2];
rz(-1.4725497) q[2];
sx q[2];
rz(2.8842928) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2858833) q[1];
sx q[1];
rz(-0.51953379) q[1];
sx q[1];
rz(-0.31182162) q[1];
x q[2];
rz(0.51185913) q[3];
sx q[3];
rz(-2.8843237) q[3];
sx q[3];
rz(-0.23990897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8740497) q[2];
sx q[2];
rz(-1.1288613) q[2];
sx q[2];
rz(1.6740602) q[2];
rz(-0.90557805) q[3];
sx q[3];
rz(-1.0970683) q[3];
sx q[3];
rz(-2.9300698) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31731376) q[0];
sx q[0];
rz(-2.2792094) q[0];
sx q[0];
rz(-0.30340075) q[0];
rz(-0.41269451) q[1];
sx q[1];
rz(-1.7957325) q[1];
sx q[1];
rz(-0.59248286) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45771407) q[0];
sx q[0];
rz(-1.787583) q[0];
sx q[0];
rz(0.44332645) q[0];
rz(1.3923774) q[2];
sx q[2];
rz(-2.1270299) q[2];
sx q[2];
rz(-1.8700333) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.35960782) q[1];
sx q[1];
rz(-2.4949412) q[1];
sx q[1];
rz(2.9293865) q[1];
x q[2];
rz(-1.1275702) q[3];
sx q[3];
rz(-0.89689287) q[3];
sx q[3];
rz(1.8473089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6233643) q[2];
sx q[2];
rz(-2.6609504) q[2];
sx q[2];
rz(0.19409689) q[2];
rz(-2.569681) q[3];
sx q[3];
rz(-1.5827936) q[3];
sx q[3];
rz(1.5552049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4647575) q[0];
sx q[0];
rz(-1.06523) q[0];
sx q[0];
rz(-1.42365) q[0];
rz(0.33547297) q[1];
sx q[1];
rz(-2.5392541) q[1];
sx q[1];
rz(-2.4981892) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5945054) q[0];
sx q[0];
rz(-0.79883655) q[0];
sx q[0];
rz(0.79669768) q[0];
rz(-1.6646447) q[2];
sx q[2];
rz(-0.99798087) q[2];
sx q[2];
rz(1.4617355) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5646518) q[1];
sx q[1];
rz(-1.1168861) q[1];
sx q[1];
rz(-2.3481525) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0213177) q[3];
sx q[3];
rz(-1.4304597) q[3];
sx q[3];
rz(-0.72718638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.246835) q[2];
sx q[2];
rz(-2.0474696) q[2];
sx q[2];
rz(-0.82478729) q[2];
rz(0.60683933) q[3];
sx q[3];
rz(-1.9347128) q[3];
sx q[3];
rz(-1.3185893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7428699) q[0];
sx q[0];
rz(-2.5720808) q[0];
sx q[0];
rz(-2.8745765) q[0];
rz(-2.6737402) q[1];
sx q[1];
rz(-2.0303969) q[1];
sx q[1];
rz(-1.9652479) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6716135) q[0];
sx q[0];
rz(-1.592127) q[0];
sx q[0];
rz(-1.6213378) q[0];
rz(0.89428561) q[2];
sx q[2];
rz(-0.91716847) q[2];
sx q[2];
rz(-2.393347) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1770498) q[1];
sx q[1];
rz(-1.6262883) q[1];
sx q[1];
rz(-2.0216136) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1517046) q[3];
sx q[3];
rz(-0.93686371) q[3];
sx q[3];
rz(-2.9798199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10845575) q[2];
sx q[2];
rz(-1.5507853) q[2];
sx q[2];
rz(0.94672686) q[2];
rz(1.4416134) q[3];
sx q[3];
rz(-1.0333034) q[3];
sx q[3];
rz(1.7472964) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93447584) q[0];
sx q[0];
rz(-1.8659135) q[0];
sx q[0];
rz(1.3502655) q[0];
rz(-2.4840202) q[1];
sx q[1];
rz(-1.5788199) q[1];
sx q[1];
rz(-1.2616166) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63245693) q[0];
sx q[0];
rz(-0.92210356) q[0];
sx q[0];
rz(-0.52059569) q[0];
x q[1];
rz(0.71102755) q[2];
sx q[2];
rz(-0.8128574) q[2];
sx q[2];
rz(-2.6136398) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4289866) q[1];
sx q[1];
rz(-1.3416051) q[1];
sx q[1];
rz(-1.5043753) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7827949) q[3];
sx q[3];
rz(-1.3984716) q[3];
sx q[3];
rz(-1.5197157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51800805) q[2];
sx q[2];
rz(-1.1058747) q[2];
sx q[2];
rz(-2.8688431) q[2];
rz(0.92877156) q[3];
sx q[3];
rz(-0.85799587) q[3];
sx q[3];
rz(-1.6509008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6269094) q[0];
sx q[0];
rz(-1.0310443) q[0];
sx q[0];
rz(0.92460257) q[0];
rz(-0.54479105) q[1];
sx q[1];
rz(-1.3489312) q[1];
sx q[1];
rz(0.10442385) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2685933) q[0];
sx q[0];
rz(-1.6571123) q[0];
sx q[0];
rz(-1.6593133) q[0];
rz(2.7265276) q[2];
sx q[2];
rz(-1.7156148) q[2];
sx q[2];
rz(0.96127779) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1510091) q[1];
sx q[1];
rz(-1.2663453) q[1];
sx q[1];
rz(-0.23465381) q[1];
x q[2];
rz(-0.35108836) q[3];
sx q[3];
rz(-1.9410053) q[3];
sx q[3];
rz(0.58634392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.32120785) q[2];
sx q[2];
rz(-1.2106004) q[2];
sx q[2];
rz(1.661181) q[2];
rz(0.004247578) q[3];
sx q[3];
rz(-2.28528) q[3];
sx q[3];
rz(2.492304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55411196) q[0];
sx q[0];
rz(-3.0453747) q[0];
sx q[0];
rz(-1.2531248) q[0];
rz(2.9907277) q[1];
sx q[1];
rz(-2.3020709) q[1];
sx q[1];
rz(1.8419267) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74053451) q[0];
sx q[0];
rz(-2.5098233) q[0];
sx q[0];
rz(0.44332544) q[0];
rz(-pi) q[1];
rz(2.0093654) q[2];
sx q[2];
rz(-1.4231496) q[2];
sx q[2];
rz(-2.4327337) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7598301) q[1];
sx q[1];
rz(-2.2564133) q[1];
sx q[1];
rz(1.5515224) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9261041) q[3];
sx q[3];
rz(-1.7607712) q[3];
sx q[3];
rz(0.097502891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8396987) q[2];
sx q[2];
rz(-1.3816741) q[2];
sx q[2];
rz(-1.985644) q[2];
rz(2.4902952) q[3];
sx q[3];
rz(-0.39172253) q[3];
sx q[3];
rz(2.7294066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10525178) q[0];
sx q[0];
rz(-1.3407433) q[0];
sx q[0];
rz(2.3403008) q[0];
rz(-0.85805145) q[1];
sx q[1];
rz(-2.4627204) q[1];
sx q[1];
rz(-0.37989315) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1042874) q[0];
sx q[0];
rz(-1.3905449) q[0];
sx q[0];
rz(2.0478115) q[0];
rz(-2.1667266) q[2];
sx q[2];
rz(-2.9306378) q[2];
sx q[2];
rz(2.4663004) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97723865) q[1];
sx q[1];
rz(-1.2697446) q[1];
sx q[1];
rz(-0.85268271) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5809459) q[3];
sx q[3];
rz(-1.1296318) q[3];
sx q[3];
rz(3.0574034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.51565591) q[2];
sx q[2];
rz(-2.8069324) q[2];
sx q[2];
rz(-2.3626309) q[2];
rz(-0.37911478) q[3];
sx q[3];
rz(-1.4232676) q[3];
sx q[3];
rz(-3.0558705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15445408) q[0];
sx q[0];
rz(-1.7310646) q[0];
sx q[0];
rz(-2.6692303) q[0];
rz(0.27913276) q[1];
sx q[1];
rz(-1.0529073) q[1];
sx q[1];
rz(2.7739024) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3084891) q[0];
sx q[0];
rz(-1.27854) q[0];
sx q[0];
rz(-1.9394642) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6515171) q[2];
sx q[2];
rz(-2.6948409) q[2];
sx q[2];
rz(1.2778145) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.93802261) q[1];
sx q[1];
rz(-1.9782776) q[1];
sx q[1];
rz(-1.8304871) q[1];
rz(-2.3261072) q[3];
sx q[3];
rz(-2.1187821) q[3];
sx q[3];
rz(-0.54875916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8184066) q[2];
sx q[2];
rz(-2.0995993) q[2];
sx q[2];
rz(2.6673178) q[2];
rz(-2.2222774) q[3];
sx q[3];
rz(-1.5654469) q[3];
sx q[3];
rz(1.7491755) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11488386) q[0];
sx q[0];
rz(-1.8395431) q[0];
sx q[0];
rz(-1.607847) q[0];
rz(-2.907091) q[1];
sx q[1];
rz(-1.0954183) q[1];
sx q[1];
rz(-0.29874994) q[1];
rz(-1.4316097) q[2];
sx q[2];
rz(-0.57568204) q[2];
sx q[2];
rz(-0.26402146) q[2];
rz(0.017757105) q[3];
sx q[3];
rz(-1.6203625) q[3];
sx q[3];
rz(0.41858471) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
