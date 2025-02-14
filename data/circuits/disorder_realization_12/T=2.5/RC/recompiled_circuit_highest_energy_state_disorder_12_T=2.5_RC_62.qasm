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
rz(-2.7589594) q[1];
sx q[1];
rz(-1.1525947) q[1];
sx q[1];
rz(-1.1421854) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4948311) q[0];
sx q[0];
rz(-1.1755623) q[0];
sx q[0];
rz(-3.0656205) q[0];
rz(-pi) q[1];
rz(-0.76337645) q[2];
sx q[2];
rz(-2.4912611) q[2];
sx q[2];
rz(1.2987607) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9153629) q[1];
sx q[1];
rz(-0.94704506) q[1];
sx q[1];
rz(-0.28304328) q[1];
rz(-0.19742404) q[3];
sx q[3];
rz(-0.89251489) q[3];
sx q[3];
rz(-1.2048282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4932058) q[2];
sx q[2];
rz(-1.1492665) q[2];
sx q[2];
rz(-0.0041848103) q[2];
rz(0.74797136) q[3];
sx q[3];
rz(-2.3145521) q[3];
sx q[3];
rz(-0.23676693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0325322) q[0];
sx q[0];
rz(-2.2219658) q[0];
sx q[0];
rz(3.076886) q[0];
rz(2.0789355) q[1];
sx q[1];
rz(-2.2149142) q[1];
sx q[1];
rz(1.0437171) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017634246) q[0];
sx q[0];
rz(-1.7746266) q[0];
sx q[0];
rz(-0.47306008) q[0];
rz(-pi) q[1];
rz(1.8383726) q[2];
sx q[2];
rz(-0.10188411) q[2];
sx q[2];
rz(0.52620483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8557094) q[1];
sx q[1];
rz(-0.51953379) q[1];
sx q[1];
rz(0.31182162) q[1];
rz(0.51185913) q[3];
sx q[3];
rz(-0.25726899) q[3];
sx q[3];
rz(0.23990897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.26754293) q[2];
sx q[2];
rz(-2.0127313) q[2];
sx q[2];
rz(-1.4675325) q[2];
rz(-2.2360146) q[3];
sx q[3];
rz(-2.0445243) q[3];
sx q[3];
rz(-2.9300698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8242789) q[0];
sx q[0];
rz(-0.86238328) q[0];
sx q[0];
rz(0.30340075) q[0];
rz(-2.7288981) q[1];
sx q[1];
rz(-1.3458601) q[1];
sx q[1];
rz(-0.59248286) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45771407) q[0];
sx q[0];
rz(-1.3540097) q[0];
sx q[0];
rz(-2.6982662) q[0];
rz(-2.5781588) q[2];
sx q[2];
rz(-1.4195002) q[2];
sx q[2];
rz(-2.7474257) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5183181) q[1];
sx q[1];
rz(-2.2006196) q[1];
sx q[1];
rz(1.4131143) q[1];
rz(-pi) q[2];
rz(2.0140225) q[3];
sx q[3];
rz(-0.89689287) q[3];
sx q[3];
rz(-1.2942838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6233643) q[2];
sx q[2];
rz(-2.6609504) q[2];
sx q[2];
rz(-2.9474958) q[2];
rz(0.57191166) q[3];
sx q[3];
rz(-1.5827936) q[3];
sx q[3];
rz(-1.5863878) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67683515) q[0];
sx q[0];
rz(-2.0763626) q[0];
sx q[0];
rz(1.42365) q[0];
rz(-2.8061197) q[1];
sx q[1];
rz(-2.5392541) q[1];
sx q[1];
rz(-2.4981892) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5945054) q[0];
sx q[0];
rz(-2.3427561) q[0];
sx q[0];
rz(0.79669768) q[0];
x q[1];
rz(2.9973028) q[2];
sx q[2];
rz(-0.57960287) q[2];
sx q[2];
rz(1.50791) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5646518) q[1];
sx q[1];
rz(-1.1168861) q[1];
sx q[1];
rz(-2.3481525) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8845064) q[3];
sx q[3];
rz(-0.47043741) q[3];
sx q[3];
rz(1.1251118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.246835) q[2];
sx q[2];
rz(-2.0474696) q[2];
sx q[2];
rz(-2.3168054) q[2];
rz(-2.5347533) q[3];
sx q[3];
rz(-1.9347128) q[3];
sx q[3];
rz(-1.3185893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3987228) q[0];
sx q[0];
rz(-2.5720808) q[0];
sx q[0];
rz(-2.8745765) q[0];
rz(-0.46785242) q[1];
sx q[1];
rz(-2.0303969) q[1];
sx q[1];
rz(1.9652479) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29821324) q[0];
sx q[0];
rz(-0.054854782) q[0];
sx q[0];
rz(1.1712267) q[0];
rz(-pi) q[1];
rz(2.3651334) q[2];
sx q[2];
rz(-1.0506223) q[2];
sx q[2];
rz(0.3683683) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96454285) q[1];
sx q[1];
rz(-1.5153043) q[1];
sx q[1];
rz(-1.1199791) q[1];
x q[2];
rz(2.1517046) q[3];
sx q[3];
rz(-2.2047289) q[3];
sx q[3];
rz(2.9798199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0331369) q[2];
sx q[2];
rz(-1.5908073) q[2];
sx q[2];
rz(2.1948658) q[2];
rz(1.6999792) q[3];
sx q[3];
rz(-1.0333034) q[3];
sx q[3];
rz(-1.7472964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(2.2071168) q[0];
sx q[0];
rz(-1.2756791) q[0];
sx q[0];
rz(-1.3502655) q[0];
rz(0.65757242) q[1];
sx q[1];
rz(-1.5627728) q[1];
sx q[1];
rz(1.2616166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63245693) q[0];
sx q[0];
rz(-0.92210356) q[0];
sx q[0];
rz(2.620997) q[0];
x q[1];
rz(-2.4305651) q[2];
sx q[2];
rz(-0.8128574) q[2];
sx q[2];
rz(-2.6136398) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2985136) q[1];
sx q[1];
rz(-1.5061146) q[1];
sx q[1];
rz(-2.9119125) q[1];
rz(-pi) q[2];
rz(-2.9654019) q[3];
sx q[3];
rz(-1.7796081) q[3];
sx q[3];
rz(0.01419078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51800805) q[2];
sx q[2];
rz(-2.035718) q[2];
sx q[2];
rz(0.27274954) q[2];
rz(0.92877156) q[3];
sx q[3];
rz(-0.85799587) q[3];
sx q[3];
rz(-1.6509008) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6269094) q[0];
sx q[0];
rz(-1.0310443) q[0];
sx q[0];
rz(-2.2169901) q[0];
rz(2.5968016) q[1];
sx q[1];
rz(-1.7926615) q[1];
sx q[1];
rz(-0.10442385) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4686779) q[0];
sx q[0];
rz(-0.12355655) q[0];
sx q[0];
rz(0.79609032) q[0];
rz(0.41506501) q[2];
sx q[2];
rz(-1.4259778) q[2];
sx q[2];
rz(0.96127779) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.663781) q[1];
sx q[1];
rz(-0.38214499) q[1];
sx q[1];
rz(-2.2078329) q[1];
rz(0.35108836) q[3];
sx q[3];
rz(-1.9410053) q[3];
sx q[3];
rz(2.5552487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8203848) q[2];
sx q[2];
rz(-1.9309923) q[2];
sx q[2];
rz(1.4804117) q[2];
rz(3.1373451) q[3];
sx q[3];
rz(-2.28528) q[3];
sx q[3];
rz(0.64928865) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5874807) q[0];
sx q[0];
rz(-3.0453747) q[0];
sx q[0];
rz(-1.8884678) q[0];
rz(2.9907277) q[1];
sx q[1];
rz(-0.8395218) q[1];
sx q[1];
rz(1.2996659) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74053451) q[0];
sx q[0];
rz(-2.5098233) q[0];
sx q[0];
rz(0.44332544) q[0];
rz(-1.9076882) q[2];
sx q[2];
rz(-0.46122069) q[2];
sx q[2];
rz(-2.5835844) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35132699) q[1];
sx q[1];
rz(-2.4557487) q[1];
sx q[1];
rz(-3.1180361) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9392936) q[3];
sx q[3];
rz(-1.222156) q[3];
sx q[3];
rz(1.5983456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.8396987) q[2];
sx q[2];
rz(-1.7599186) q[2];
sx q[2];
rz(1.1559486) q[2];
rz(2.4902952) q[3];
sx q[3];
rz(-0.39172253) q[3];
sx q[3];
rz(2.7294066) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0363409) q[0];
sx q[0];
rz(-1.3407433) q[0];
sx q[0];
rz(-0.80129188) q[0];
rz(-0.85805145) q[1];
sx q[1];
rz(-2.4627204) q[1];
sx q[1];
rz(2.7616995) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86738619) q[0];
sx q[0];
rz(-0.50747061) q[0];
sx q[0];
rz(1.9486289) q[0];
rz(3.0219737) q[2];
sx q[2];
rz(-1.3966171) q[2];
sx q[2];
rz(1.8599138) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2945248) q[1];
sx q[1];
rz(-2.2502568) q[1];
sx q[1];
rz(-2.7505198) q[1];
x q[2];
rz(1.5809459) q[3];
sx q[3];
rz(-1.1296318) q[3];
sx q[3];
rz(3.0574034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51565591) q[2];
sx q[2];
rz(-0.3346602) q[2];
sx q[2];
rz(-0.77896172) q[2];
rz(2.7624779) q[3];
sx q[3];
rz(-1.4232676) q[3];
sx q[3];
rz(0.085722119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15445408) q[0];
sx q[0];
rz(-1.4105281) q[0];
sx q[0];
rz(2.6692303) q[0];
rz(-0.27913276) q[1];
sx q[1];
rz(-1.0529073) q[1];
sx q[1];
rz(-2.7739024) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5147509) q[0];
sx q[0];
rz(-1.923133) q[0];
sx q[0];
rz(2.8295838) q[0];
x q[1];
rz(1.3490178) q[2];
sx q[2];
rz(-1.1797172) q[2];
sx q[2];
rz(1.329601) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.73767978) q[1];
sx q[1];
rz(-1.8087937) q[1];
sx q[1];
rz(0.42003553) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81548549) q[3];
sx q[3];
rz(-2.1187821) q[3];
sx q[3];
rz(2.5928335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8184066) q[2];
sx q[2];
rz(-1.0419934) q[2];
sx q[2];
rz(-2.6673178) q[2];
rz(0.91931528) q[3];
sx q[3];
rz(-1.5761458) q[3];
sx q[3];
rz(-1.7491755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11488386) q[0];
sx q[0];
rz(-1.8395431) q[0];
sx q[0];
rz(-1.607847) q[0];
rz(0.23450163) q[1];
sx q[1];
rz(-1.0954183) q[1];
sx q[1];
rz(-0.29874994) q[1];
rz(2.1420494) q[2];
sx q[2];
rz(-1.4951946) q[2];
sx q[2];
rz(-1.9517938) q[2];
rz(1.5212223) q[3];
sx q[3];
rz(-1.5885316) q[3];
sx q[3];
rz(1.9902609) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
