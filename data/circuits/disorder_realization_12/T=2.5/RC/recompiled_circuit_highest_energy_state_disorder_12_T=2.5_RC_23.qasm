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
rz(-1.9528376) q[0];
sx q[0];
rz(-0.83851695) q[0];
sx q[0];
rz(-1.4128348) q[0];
rz(0.38263327) q[1];
sx q[1];
rz(-1.9889979) q[1];
sx q[1];
rz(1.1421854) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6467616) q[0];
sx q[0];
rz(-1.9660304) q[0];
sx q[0];
rz(3.0656205) q[0];
rz(-0.5025592) q[2];
sx q[2];
rz(-1.1389073) q[2];
sx q[2];
rz(-2.2186861) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.22622977) q[1];
sx q[1];
rz(-0.94704506) q[1];
sx q[1];
rz(-0.28304328) q[1];
x q[2];
rz(-1.3320311) q[3];
sx q[3];
rz(-2.4395535) q[3];
sx q[3];
rz(-2.2453737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6483868) q[2];
sx q[2];
rz(-1.1492665) q[2];
sx q[2];
rz(0.0041848103) q[2];
rz(-0.74797136) q[3];
sx q[3];
rz(-0.82704058) q[3];
sx q[3];
rz(2.9048257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0325322) q[0];
sx q[0];
rz(-0.91962686) q[0];
sx q[0];
rz(3.076886) q[0];
rz(2.0789355) q[1];
sx q[1];
rz(-2.2149142) q[1];
sx q[1];
rz(-2.0978755) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4499265) q[0];
sx q[0];
rz(-2.0332893) q[0];
sx q[0];
rz(1.7989547) q[0];
rz(-pi) q[1];
rz(1.4725141) q[2];
sx q[2];
rz(-1.5976906) q[2];
sx q[2];
rz(-1.310845) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5001561) q[1];
sx q[1];
rz(-2.0629597) q[1];
sx q[1];
rz(-1.3970966) q[1];
rz(-pi) q[2];
rz(0.51185913) q[3];
sx q[3];
rz(-2.8843237) q[3];
sx q[3];
rz(2.9016837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8740497) q[2];
sx q[2];
rz(-2.0127313) q[2];
sx q[2];
rz(-1.6740602) q[2];
rz(2.2360146) q[3];
sx q[3];
rz(-2.0445243) q[3];
sx q[3];
rz(-0.21152285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8242789) q[0];
sx q[0];
rz(-2.2792094) q[0];
sx q[0];
rz(2.8381919) q[0];
rz(-0.41269451) q[1];
sx q[1];
rz(-1.3458601) q[1];
sx q[1];
rz(0.59248286) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68776206) q[0];
sx q[0];
rz(-0.49031165) q[0];
sx q[0];
rz(-2.6672336) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3923774) q[2];
sx q[2];
rz(-1.0145628) q[2];
sx q[2];
rz(-1.2715593) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.62327451) q[1];
sx q[1];
rz(-2.2006196) q[1];
sx q[1];
rz(-1.7284784) q[1];
x q[2];
rz(0.49280102) q[3];
sx q[3];
rz(-2.3544724) q[3];
sx q[3];
rz(-1.9446179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6233643) q[2];
sx q[2];
rz(-2.6609504) q[2];
sx q[2];
rz(2.9474958) q[2];
rz(0.57191166) q[3];
sx q[3];
rz(-1.5587991) q[3];
sx q[3];
rz(-1.5552049) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4647575) q[0];
sx q[0];
rz(-2.0763626) q[0];
sx q[0];
rz(1.7179426) q[0];
rz(-0.33547297) q[1];
sx q[1];
rz(-2.5392541) q[1];
sx q[1];
rz(-0.64340341) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62221974) q[0];
sx q[0];
rz(-1.046139) q[0];
sx q[0];
rz(2.2043246) q[0];
rz(1.6646447) q[2];
sx q[2];
rz(-0.99798087) q[2];
sx q[2];
rz(1.6798572) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.58650342) q[1];
sx q[1];
rz(-0.88876206) q[1];
sx q[1];
rz(-2.5413496) q[1];
x q[2];
rz(-2.9859391) q[3];
sx q[3];
rz(-2.0165679) q[3];
sx q[3];
rz(0.77605493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.246835) q[2];
sx q[2];
rz(-2.0474696) q[2];
sx q[2];
rz(-2.3168054) q[2];
rz(0.60683933) q[3];
sx q[3];
rz(-1.2068799) q[3];
sx q[3];
rz(1.3185893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7428699) q[0];
sx q[0];
rz(-2.5720808) q[0];
sx q[0];
rz(2.8745765) q[0];
rz(0.46785242) q[1];
sx q[1];
rz(-1.1111958) q[1];
sx q[1];
rz(-1.1763447) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1018961) q[0];
sx q[0];
rz(-1.6213263) q[0];
sx q[0];
rz(-3.1202348) q[0];
rz(-pi) q[1];
rz(2.247307) q[2];
sx q[2];
rz(-2.2244242) q[2];
sx q[2];
rz(0.7482457) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5084978) q[1];
sx q[1];
rz(-2.0208686) q[1];
sx q[1];
rz(-0.061636713) q[1];
x q[2];
rz(0.98988804) q[3];
sx q[3];
rz(-0.93686371) q[3];
sx q[3];
rz(2.9798199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0331369) q[2];
sx q[2];
rz(-1.5507853) q[2];
sx q[2];
rz(-0.94672686) q[2];
rz(-1.4416134) q[3];
sx q[3];
rz(-1.0333034) q[3];
sx q[3];
rz(1.3942963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93447584) q[0];
sx q[0];
rz(-1.2756791) q[0];
sx q[0];
rz(-1.7913272) q[0];
rz(0.65757242) q[1];
sx q[1];
rz(-1.5627728) q[1];
sx q[1];
rz(1.2616166) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3917106) q[0];
sx q[0];
rz(-0.80750033) q[0];
sx q[0];
rz(-2.1514284) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4665537) q[2];
sx q[2];
rz(-1.0770105) q[2];
sx q[2];
rz(1.5774278) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4289866) q[1];
sx q[1];
rz(-1.7999876) q[1];
sx q[1];
rz(1.5043753) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17619074) q[3];
sx q[3];
rz(-1.3619845) q[3];
sx q[3];
rz(0.01419078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51800805) q[2];
sx q[2];
rz(-1.1058747) q[2];
sx q[2];
rz(-0.27274954) q[2];
rz(0.92877156) q[3];
sx q[3];
rz(-0.85799587) q[3];
sx q[3];
rz(-1.6509008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.6269094) q[0];
sx q[0];
rz(-2.1105483) q[0];
sx q[0];
rz(-2.2169901) q[0];
rz(0.54479105) q[1];
sx q[1];
rz(-1.7926615) q[1];
sx q[1];
rz(0.10442385) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4514465) q[0];
sx q[0];
rz(-1.4826097) q[0];
sx q[0];
rz(3.0549391) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7265276) q[2];
sx q[2];
rz(-1.4259778) q[2];
sx q[2];
rz(-0.96127779) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.663781) q[1];
sx q[1];
rz(-0.38214499) q[1];
sx q[1];
rz(2.2078329) q[1];
x q[2];
rz(-1.9627294) q[3];
sx q[3];
rz(-1.8971803) q[3];
sx q[3];
rz(1.1162023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.32120785) q[2];
sx q[2];
rz(-1.2106004) q[2];
sx q[2];
rz(-1.4804117) q[2];
rz(-0.004247578) q[3];
sx q[3];
rz(-2.28528) q[3];
sx q[3];
rz(0.64928865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55411196) q[0];
sx q[0];
rz(-0.09621796) q[0];
sx q[0];
rz(-1.2531248) q[0];
rz(2.9907277) q[1];
sx q[1];
rz(-2.3020709) q[1];
sx q[1];
rz(-1.2996659) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2723994) q[0];
sx q[0];
rz(-2.1335094) q[0];
sx q[0];
rz(-1.8749694) q[0];
rz(-0.1628218) q[2];
sx q[2];
rz(-1.1373242) q[2];
sx q[2];
rz(-0.93082817) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.35132699) q[1];
sx q[1];
rz(-0.68584397) q[1];
sx q[1];
rz(3.1180361) q[1];
rz(-pi) q[2];
rz(-1.2154886) q[3];
sx q[3];
rz(-1.3808215) q[3];
sx q[3];
rz(0.097502891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8396987) q[2];
sx q[2];
rz(-1.7599186) q[2];
sx q[2];
rz(1.1559486) q[2];
rz(-0.65129748) q[3];
sx q[3];
rz(-0.39172253) q[3];
sx q[3];
rz(-0.41218606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0363409) q[0];
sx q[0];
rz(-1.3407433) q[0];
sx q[0];
rz(-2.3403008) q[0];
rz(-0.85805145) q[1];
sx q[1];
rz(-0.67887226) q[1];
sx q[1];
rz(-2.7616995) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1042874) q[0];
sx q[0];
rz(-1.7510478) q[0];
sx q[0];
rz(-2.0478115) q[0];
rz(0.97486603) q[2];
sx q[2];
rz(-2.9306378) q[2];
sx q[2];
rz(2.4663004) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8752138) q[1];
sx q[1];
rz(-2.3734178) q[1];
sx q[1];
rz(-1.1298864) q[1];
rz(-1.5809459) q[3];
sx q[3];
rz(-2.0119609) q[3];
sx q[3];
rz(3.0574034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51565591) q[2];
sx q[2];
rz(-2.8069324) q[2];
sx q[2];
rz(-0.77896172) q[2];
rz(-0.37911478) q[3];
sx q[3];
rz(-1.4232676) q[3];
sx q[3];
rz(-3.0558705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15445408) q[0];
sx q[0];
rz(-1.4105281) q[0];
sx q[0];
rz(-0.4723624) q[0];
rz(0.27913276) q[1];
sx q[1];
rz(-1.0529073) q[1];
sx q[1];
rz(-0.36769029) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5147509) q[0];
sx q[0];
rz(-1.923133) q[0];
sx q[0];
rz(-2.8295838) q[0];
rz(-1.7925749) q[2];
sx q[2];
rz(-1.9618755) q[2];
sx q[2];
rz(1.8119916) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7941741) q[1];
sx q[1];
rz(-2.6623335) q[1];
sx q[1];
rz(-0.53665953) q[1];
rz(-0.69766694) q[3];
sx q[3];
rz(-2.1958477) q[3];
sx q[3];
rz(-0.56599697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8184066) q[2];
sx q[2];
rz(-1.0419934) q[2];
sx q[2];
rz(-0.47427487) q[2];
rz(2.2222774) q[3];
sx q[3];
rz(-1.5761458) q[3];
sx q[3];
rz(-1.3924172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0267088) q[0];
sx q[0];
rz(-1.8395431) q[0];
sx q[0];
rz(-1.607847) q[0];
rz(2.907091) q[1];
sx q[1];
rz(-2.0461743) q[1];
sx q[1];
rz(2.8428427) q[1];
rz(-3.0517921) q[2];
sx q[2];
rz(-2.1402146) q[2];
sx q[2];
rz(-0.42950523) q[2];
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
