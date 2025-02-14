OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.65052819) q[0];
sx q[0];
rz(-1.0125546) q[0];
sx q[0];
rz(0.92232409) q[0];
rz(-0.84696472) q[1];
sx q[1];
rz(-1.6672517) q[1];
sx q[1];
rz(0.21811952) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51061741) q[0];
sx q[0];
rz(-1.6625064) q[0];
sx q[0];
rz(2.6789078) q[0];
x q[1];
rz(0.37801554) q[2];
sx q[2];
rz(-1.6089919) q[2];
sx q[2];
rz(0.60249828) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0099758) q[1];
sx q[1];
rz(-0.22356114) q[1];
sx q[1];
rz(1.7625336) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7401198) q[3];
sx q[3];
rz(-0.43531951) q[3];
sx q[3];
rz(0.12745133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.073079022) q[2];
sx q[2];
rz(-1.2129236) q[2];
sx q[2];
rz(-0.53654137) q[2];
rz(1.0088629) q[3];
sx q[3];
rz(-2.3705685) q[3];
sx q[3];
rz(2.3276276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5040078) q[0];
sx q[0];
rz(-1.8300087) q[0];
sx q[0];
rz(3.1245533) q[0];
rz(-1.0631961) q[1];
sx q[1];
rz(-2.3737962) q[1];
sx q[1];
rz(2.1868736) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9794036) q[0];
sx q[0];
rz(-1.5155751) q[0];
sx q[0];
rz(2.9121141) q[0];
x q[1];
rz(-0.164523) q[2];
sx q[2];
rz(-1.8662819) q[2];
sx q[2];
rz(0.72876677) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0481457) q[1];
sx q[1];
rz(-2.8110782) q[1];
sx q[1];
rz(2.3814209) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.019308643) q[3];
sx q[3];
rz(-1.5175765) q[3];
sx q[3];
rz(1.9198708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.22975989) q[2];
sx q[2];
rz(-2.2288897) q[2];
sx q[2];
rz(2.0274053) q[2];
rz(-1.1344502) q[3];
sx q[3];
rz(-0.19616923) q[3];
sx q[3];
rz(-2.9451059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5512307) q[0];
sx q[0];
rz(-2.1901665) q[0];
sx q[0];
rz(-0.1828585) q[0];
rz(0.081347801) q[1];
sx q[1];
rz(-2.3902939) q[1];
sx q[1];
rz(-2.523211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9524117) q[0];
sx q[0];
rz(-0.78972406) q[0];
sx q[0];
rz(0.10297063) q[0];
x q[1];
rz(0.19962991) q[2];
sx q[2];
rz(-2.3702894) q[2];
sx q[2];
rz(-0.16801258) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8152541) q[1];
sx q[1];
rz(-0.36838057) q[1];
sx q[1];
rz(-1.9659564) q[1];
x q[2];
rz(0.62348714) q[3];
sx q[3];
rz(-2.1716431) q[3];
sx q[3];
rz(-0.86756016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0028093) q[2];
sx q[2];
rz(-1.8178136) q[2];
sx q[2];
rz(-2.7749824) q[2];
rz(-1.9256516) q[3];
sx q[3];
rz(-1.1953019) q[3];
sx q[3];
rz(2.478157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4902896) q[0];
sx q[0];
rz(-1.9090575) q[0];
sx q[0];
rz(-2.41462) q[0];
rz(0.90826774) q[1];
sx q[1];
rz(-0.5779225) q[1];
sx q[1];
rz(2.0982826) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9350727) q[0];
sx q[0];
rz(-1.5653725) q[0];
sx q[0];
rz(2.4804153) q[0];
rz(-pi) q[1];
rz(-1.415148) q[2];
sx q[2];
rz(-2.7619432) q[2];
sx q[2];
rz(1.8774892) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2542782) q[1];
sx q[1];
rz(-1.320854) q[1];
sx q[1];
rz(-1.1123841) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0049694) q[3];
sx q[3];
rz(-2.5006066) q[3];
sx q[3];
rz(2.084465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4116481) q[2];
sx q[2];
rz(-2.8204155) q[2];
sx q[2];
rz(1.8357065) q[2];
rz(0.11792396) q[3];
sx q[3];
rz(-0.4685466) q[3];
sx q[3];
rz(-1.0606631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0333198) q[0];
sx q[0];
rz(-2.1554027) q[0];
sx q[0];
rz(1.6075851) q[0];
rz(0.29574212) q[1];
sx q[1];
rz(-1.5402126) q[1];
sx q[1];
rz(1.6827513) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5285501) q[0];
sx q[0];
rz(-0.82697059) q[0];
sx q[0];
rz(3.0463329) q[0];
x q[1];
rz(-1.940002) q[2];
sx q[2];
rz(-1.8593238) q[2];
sx q[2];
rz(2.8099589) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70965687) q[1];
sx q[1];
rz(-1.767358) q[1];
sx q[1];
rz(-2.2811969) q[1];
x q[2];
rz(-2.1516861) q[3];
sx q[3];
rz(-1.3767929) q[3];
sx q[3];
rz(-1.7542183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7532588) q[2];
sx q[2];
rz(-0.39187852) q[2];
sx q[2];
rz(-2.4957116) q[2];
rz(-1.9258457) q[3];
sx q[3];
rz(-1.4589717) q[3];
sx q[3];
rz(-1.7997883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.129824) q[0];
sx q[0];
rz(-2.3079066) q[0];
sx q[0];
rz(-0.60761333) q[0];
rz(1.1786849) q[1];
sx q[1];
rz(-1.8557529) q[1];
sx q[1];
rz(2.7778621) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0775438) q[0];
sx q[0];
rz(-2.0850967) q[0];
sx q[0];
rz(-0.27405996) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9035788) q[2];
sx q[2];
rz(-1.5269035) q[2];
sx q[2];
rz(1.8157168) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1235001) q[1];
sx q[1];
rz(-0.82595968) q[1];
sx q[1];
rz(-2.9781746) q[1];
rz(-pi) q[2];
rz(-3.1167459) q[3];
sx q[3];
rz(-2.1669706) q[3];
sx q[3];
rz(2.5667584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58934775) q[2];
sx q[2];
rz(-0.10422464) q[2];
sx q[2];
rz(1.1927401) q[2];
rz(-0.73181152) q[3];
sx q[3];
rz(-1.7164427) q[3];
sx q[3];
rz(1.653479) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8530387) q[0];
sx q[0];
rz(-2.9712501) q[0];
sx q[0];
rz(-1.6188251) q[0];
rz(-1.4672) q[1];
sx q[1];
rz(-1.8312981) q[1];
sx q[1];
rz(0.67749643) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.753669) q[0];
sx q[0];
rz(-2.4437527) q[0];
sx q[0];
rz(1.4927255) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67393731) q[2];
sx q[2];
rz(-0.29602414) q[2];
sx q[2];
rz(-0.37002555) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4851393) q[1];
sx q[1];
rz(-2.5514126) q[1];
sx q[1];
rz(-2.3132669) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17685276) q[3];
sx q[3];
rz(-0.74027432) q[3];
sx q[3];
rz(-2.3725584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39685321) q[2];
sx q[2];
rz(-2.2262959) q[2];
sx q[2];
rz(1.8358561) q[2];
rz(0.33852494) q[3];
sx q[3];
rz(-1.1208231) q[3];
sx q[3];
rz(-0.6849851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6029538) q[0];
sx q[0];
rz(-0.21553497) q[0];
sx q[0];
rz(-0.18390528) q[0];
rz(-1.6869847) q[1];
sx q[1];
rz(-0.55211663) q[1];
sx q[1];
rz(0.49585453) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3199098) q[0];
sx q[0];
rz(-1.8368854) q[0];
sx q[0];
rz(2.0771107) q[0];
x q[1];
rz(2.3681201) q[2];
sx q[2];
rz(-2.2272553) q[2];
sx q[2];
rz(-0.17420775) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.293893) q[1];
sx q[1];
rz(-1.2155079) q[1];
sx q[1];
rz(-1.1637264) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0184278) q[3];
sx q[3];
rz(-1.6721069) q[3];
sx q[3];
rz(2.5309569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.054606525) q[2];
sx q[2];
rz(-0.84422529) q[2];
sx q[2];
rz(-1.1620713) q[2];
rz(1.453513) q[3];
sx q[3];
rz(-2.1789357) q[3];
sx q[3];
rz(1.8614205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42089713) q[0];
sx q[0];
rz(-1.9891885) q[0];
sx q[0];
rz(-0.56588093) q[0];
rz(-0.3903009) q[1];
sx q[1];
rz(-1.5696328) q[1];
sx q[1];
rz(1.8338667) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824664) q[0];
sx q[0];
rz(-2.0604366) q[0];
sx q[0];
rz(0.35207502) q[0];
rz(1.9487319) q[2];
sx q[2];
rz(-1.3728752) q[2];
sx q[2];
rz(1.6557616) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3257287) q[1];
sx q[1];
rz(-0.96885175) q[1];
sx q[1];
rz(1.25156) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89687895) q[3];
sx q[3];
rz(-2.9510636) q[3];
sx q[3];
rz(-2.6485557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2719443) q[2];
sx q[2];
rz(-2.7248236) q[2];
sx q[2];
rz(2.1040253) q[2];
rz(-1.3197673) q[3];
sx q[3];
rz(-0.82873738) q[3];
sx q[3];
rz(-0.64798361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8808402) q[0];
sx q[0];
rz(-2.1639731) q[0];
sx q[0];
rz(0.64224893) q[0];
rz(-1.8107268) q[1];
sx q[1];
rz(-2.2763177) q[1];
sx q[1];
rz(-2.9790402) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88961381) q[0];
sx q[0];
rz(-1.9601213) q[0];
sx q[0];
rz(2.2963803) q[0];
x q[1];
rz(-1.2106718) q[2];
sx q[2];
rz(-1.2136493) q[2];
sx q[2];
rz(-1.2993297) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0916718) q[1];
sx q[1];
rz(-2.0494048) q[1];
sx q[1];
rz(-2.9438627) q[1];
x q[2];
rz(0.58329669) q[3];
sx q[3];
rz(-1.6371033) q[3];
sx q[3];
rz(-0.60140677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6315397) q[2];
sx q[2];
rz(-1.2756462) q[2];
sx q[2];
rz(1.0408545) q[2];
rz(-0.96796525) q[3];
sx q[3];
rz(-2.0716045) q[3];
sx q[3];
rz(-0.66106838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80617245) q[0];
sx q[0];
rz(-1.5652884) q[0];
sx q[0];
rz(-0.096927222) q[0];
rz(0.23282911) q[1];
sx q[1];
rz(-1.0284582) q[1];
sx q[1];
rz(-3.1340541) q[1];
rz(0.54388028) q[2];
sx q[2];
rz(-0.82533045) q[2];
sx q[2];
rz(-2.5173204) q[2];
rz(1.9540174) q[3];
sx q[3];
rz(-0.90771994) q[3];
sx q[3];
rz(-0.049605443) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
