OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7856287) q[0];
sx q[0];
rz(-0.2956737) q[0];
sx q[0];
rz(9.8386968) q[0];
rz(0.26894459) q[1];
sx q[1];
rz(-0.36965951) q[1];
sx q[1];
rz(-1.2764021) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39525014) q[0];
sx q[0];
rz(-1.479404) q[0];
sx q[0];
rz(1.6618098) q[0];
rz(-pi) q[1];
rz(-1.6501313) q[2];
sx q[2];
rz(-0.66986194) q[2];
sx q[2];
rz(-1.2487433) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.947283) q[1];
sx q[1];
rz(-2.3230246) q[1];
sx q[1];
rz(1.964393) q[1];
rz(-pi) q[2];
rz(-2.5475575) q[3];
sx q[3];
rz(-1.7230936) q[3];
sx q[3];
rz(-1.9826979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.74497491) q[2];
sx q[2];
rz(-1.9048385) q[2];
sx q[2];
rz(-1.940894) q[2];
rz(0.68136224) q[3];
sx q[3];
rz(-1.6187637) q[3];
sx q[3];
rz(-0.38667005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.50614) q[0];
sx q[0];
rz(-0.30158392) q[0];
sx q[0];
rz(-0.26148456) q[0];
rz(-1.5226978) q[1];
sx q[1];
rz(-0.50917429) q[1];
sx q[1];
rz(-1.8656628) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0700713) q[0];
sx q[0];
rz(-2.1198648) q[0];
sx q[0];
rz(1.683805) q[0];
rz(-pi) q[1];
rz(1.0307113) q[2];
sx q[2];
rz(-2.1505877) q[2];
sx q[2];
rz(2.3183398) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9918242) q[1];
sx q[1];
rz(-0.64095488) q[1];
sx q[1];
rz(-2.5800406) q[1];
rz(-pi) q[2];
rz(1.719763) q[3];
sx q[3];
rz(-1.5368122) q[3];
sx q[3];
rz(-0.33992088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2168938) q[2];
sx q[2];
rz(-2.108722) q[2];
sx q[2];
rz(2.4135598) q[2];
rz(1.3701471) q[3];
sx q[3];
rz(-2.5244505) q[3];
sx q[3];
rz(-0.035577687) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034123357) q[0];
sx q[0];
rz(-1.1119482) q[0];
sx q[0];
rz(2.4566101) q[0];
rz(-2.0383539) q[1];
sx q[1];
rz(-0.61956844) q[1];
sx q[1];
rz(2.129295) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023803614) q[0];
sx q[0];
rz(-0.53142457) q[0];
sx q[0];
rz(-0.84367911) q[0];
rz(-2.2980897) q[2];
sx q[2];
rz(-0.25878227) q[2];
sx q[2];
rz(-2.7121787) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6382513) q[1];
sx q[1];
rz(-1.6972528) q[1];
sx q[1];
rz(2.9328037) q[1];
rz(-0.7638974) q[3];
sx q[3];
rz(-1.9490846) q[3];
sx q[3];
rz(-1.0920785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.181902) q[2];
sx q[2];
rz(-1.5896553) q[2];
sx q[2];
rz(2.336179) q[2];
rz(-0.81625932) q[3];
sx q[3];
rz(-0.62556848) q[3];
sx q[3];
rz(0.094154112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020298088) q[0];
sx q[0];
rz(-0.90885201) q[0];
sx q[0];
rz(-0.62776172) q[0];
rz(-3.032387) q[1];
sx q[1];
rz(-0.86995482) q[1];
sx q[1];
rz(-0.83762082) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8872466) q[0];
sx q[0];
rz(-1.6372879) q[0];
sx q[0];
rz(2.6496135) q[0];
rz(1.2455279) q[2];
sx q[2];
rz(-2.5211689) q[2];
sx q[2];
rz(-1.9577946) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7594436) q[1];
sx q[1];
rz(-1.7045705) q[1];
sx q[1];
rz(-2.4390038) q[1];
rz(-pi) q[2];
rz(-0.37552278) q[3];
sx q[3];
rz(-1.8584455) q[3];
sx q[3];
rz(-0.2660397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3179021) q[2];
sx q[2];
rz(-0.48419848) q[2];
sx q[2];
rz(2.79706) q[2];
rz(-1.7349998) q[3];
sx q[3];
rz(-0.35511261) q[3];
sx q[3];
rz(3.0936892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(2.542881) q[0];
sx q[0];
rz(-2.5051835) q[0];
sx q[0];
rz(2.6357546) q[0];
rz(2.089031) q[1];
sx q[1];
rz(-2.274175) q[1];
sx q[1];
rz(0.94598407) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1051798) q[0];
sx q[0];
rz(-0.37659058) q[0];
sx q[0];
rz(-0.74684493) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35752488) q[2];
sx q[2];
rz(-1.4357743) q[2];
sx q[2];
rz(2.6808878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.16984384) q[1];
sx q[1];
rz(-1.1976114) q[1];
sx q[1];
rz(0.62974522) q[1];
x q[2];
rz(1.9993773) q[3];
sx q[3];
rz(-0.16332291) q[3];
sx q[3];
rz(1.4455011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0537009) q[2];
sx q[2];
rz(-2.044007) q[2];
sx q[2];
rz(1.5778479) q[2];
rz(2.0569233) q[3];
sx q[3];
rz(-2.2659149) q[3];
sx q[3];
rz(-0.55115551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56187335) q[0];
sx q[0];
rz(-0.60573524) q[0];
sx q[0];
rz(-2.9015923) q[0];
rz(2.2198246) q[1];
sx q[1];
rz(-1.0488989) q[1];
sx q[1];
rz(1.1704495) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83607996) q[0];
sx q[0];
rz(-1.613998) q[0];
sx q[0];
rz(0.57010285) q[0];
rz(-pi) q[1];
rz(2.8113995) q[2];
sx q[2];
rz(-2.7666436) q[2];
sx q[2];
rz(-0.2548616) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.45615087) q[1];
sx q[1];
rz(-2.4439881) q[1];
sx q[1];
rz(0.82932034) q[1];
rz(-pi) q[2];
rz(3.0441566) q[3];
sx q[3];
rz(-1.9698922) q[3];
sx q[3];
rz(0.68828415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.74448284) q[2];
sx q[2];
rz(-2.3220671) q[2];
sx q[2];
rz(-2.5605555) q[2];
rz(0.93823141) q[3];
sx q[3];
rz(-0.56648985) q[3];
sx q[3];
rz(-2.4375622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.13643232) q[0];
sx q[0];
rz(-2.4649824) q[0];
sx q[0];
rz(-2.6410979) q[0];
rz(-2.9035134) q[1];
sx q[1];
rz(-1.6903189) q[1];
sx q[1];
rz(-2.4949825) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85486551) q[0];
sx q[0];
rz(-1.5359794) q[0];
sx q[0];
rz(-1.5690194) q[0];
x q[1];
rz(1.4562723) q[2];
sx q[2];
rz(-1.537064) q[2];
sx q[2];
rz(2.9303868) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.889786) q[1];
sx q[1];
rz(-1.6908892) q[1];
sx q[1];
rz(2.0577862) q[1];
rz(1.929333) q[3];
sx q[3];
rz(-1.6286116) q[3];
sx q[3];
rz(-1.1392347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2639192) q[2];
sx q[2];
rz(-0.71920243) q[2];
sx q[2];
rz(2.4379697) q[2];
rz(-1.8536812) q[3];
sx q[3];
rz(-2.0813007) q[3];
sx q[3];
rz(-0.50584403) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25381655) q[0];
sx q[0];
rz(-1.1340589) q[0];
sx q[0];
rz(1.3251086) q[0];
rz(-2.3690986) q[1];
sx q[1];
rz(-1.1224116) q[1];
sx q[1];
rz(0.17556369) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2529375) q[0];
sx q[0];
rz(-2.3217259) q[0];
sx q[0];
rz(-1.8273414) q[0];
x q[1];
rz(-2.2261691) q[2];
sx q[2];
rz(-0.88214126) q[2];
sx q[2];
rz(1.4530393) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0649291) q[1];
sx q[1];
rz(-0.48427072) q[1];
sx q[1];
rz(-1.8736429) q[1];
rz(-1.4453078) q[3];
sx q[3];
rz(-0.85832046) q[3];
sx q[3];
rz(-3.0355796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.10302155) q[2];
sx q[2];
rz(-2.0671637) q[2];
sx q[2];
rz(1.0938905) q[2];
rz(1.966018) q[3];
sx q[3];
rz(-2.3837377) q[3];
sx q[3];
rz(-0.28483835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4232101) q[0];
sx q[0];
rz(-3.0960313) q[0];
sx q[0];
rz(-1.3257931) q[0];
rz(-2.6153053) q[1];
sx q[1];
rz(-0.86295366) q[1];
sx q[1];
rz(-2.3525499) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5518582) q[0];
sx q[0];
rz(-1.6673256) q[0];
sx q[0];
rz(1.7328698) q[0];
rz(-pi) q[1];
rz(-0.35576917) q[2];
sx q[2];
rz(-1.3311183) q[2];
sx q[2];
rz(3.0955448) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7154768) q[1];
sx q[1];
rz(-1.9144692) q[1];
sx q[1];
rz(2.2261376) q[1];
rz(2.4715244) q[3];
sx q[3];
rz(-2.2446771) q[3];
sx q[3];
rz(1.4953227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9290756) q[2];
sx q[2];
rz(-0.93364659) q[2];
sx q[2];
rz(-0.90544256) q[2];
rz(-0.91600156) q[3];
sx q[3];
rz(-1.4429049) q[3];
sx q[3];
rz(-1.0441095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7020029) q[0];
sx q[0];
rz(-0.8332533) q[0];
sx q[0];
rz(-2.2138017) q[0];
rz(2.0480428) q[1];
sx q[1];
rz(-0.55892006) q[1];
sx q[1];
rz(3.0082026) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0693215) q[0];
sx q[0];
rz(-1.7786553) q[0];
sx q[0];
rz(-0.29671191) q[0];
x q[1];
rz(2.1422562) q[2];
sx q[2];
rz(-1.396787) q[2];
sx q[2];
rz(1.0854173) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.16001734) q[1];
sx q[1];
rz(-0.80343738) q[1];
sx q[1];
rz(2.395242) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9319256) q[3];
sx q[3];
rz(-2.0197967) q[3];
sx q[3];
rz(-2.1213139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5043958) q[2];
sx q[2];
rz(-1.7375676) q[2];
sx q[2];
rz(1.2633598) q[2];
rz(-2.0718306) q[3];
sx q[3];
rz(-2.1061149) q[3];
sx q[3];
rz(-0.076816471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1509811) q[0];
sx q[0];
rz(-0.69775109) q[0];
sx q[0];
rz(1.3569111) q[0];
rz(1.4487343) q[1];
sx q[1];
rz(-2.3285463) q[1];
sx q[1];
rz(2.9978233) q[1];
rz(-2.0781381) q[2];
sx q[2];
rz(-0.58488556) q[2];
sx q[2];
rz(1.9716138) q[2];
rz(0.95851687) q[3];
sx q[3];
rz(-2.7296441) q[3];
sx q[3];
rz(-1.7492528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
