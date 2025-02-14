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
rz(2.3866374) q[0];
sx q[0];
rz(-0.75140262) q[0];
sx q[0];
rz(-1.0487392) q[0];
rz(0.41930786) q[1];
sx q[1];
rz(-1.3860621) q[1];
sx q[1];
rz(-3.0233033) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70067507) q[0];
sx q[0];
rz(-2.5108733) q[0];
sx q[0];
rz(-1.0998902) q[0];
x q[1];
rz(-1.5871136) q[2];
sx q[2];
rz(-1.6217578) q[2];
sx q[2];
rz(-1.9555443) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.018838192) q[1];
sx q[1];
rz(-1.600482) q[1];
sx q[1];
rz(-1.5877941) q[1];
x q[2];
rz(-2.8382073) q[3];
sx q[3];
rz(-1.4974563) q[3];
sx q[3];
rz(-0.23305063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2179541) q[2];
sx q[2];
rz(-3.1380234) q[2];
sx q[2];
rz(0.16036073) q[2];
rz(-1.1637566) q[3];
sx q[3];
rz(-1.0842423) q[3];
sx q[3];
rz(-0.81419301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9775951) q[0];
sx q[0];
rz(-1.6544592) q[0];
sx q[0];
rz(2.1556222) q[0];
rz(-1.5892971) q[1];
sx q[1];
rz(-2.8542216) q[1];
sx q[1];
rz(1.5602962) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6628806) q[0];
sx q[0];
rz(-1.1887738) q[0];
sx q[0];
rz(2.0860614) q[0];
rz(-pi) q[1];
rz(-1.549717) q[2];
sx q[2];
rz(-0.97641845) q[2];
sx q[2];
rz(0.057986857) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4680921) q[1];
sx q[1];
rz(-1.7335658) q[1];
sx q[1];
rz(0.36483187) q[1];
x q[2];
rz(3.017707) q[3];
sx q[3];
rz(-2.200211) q[3];
sx q[3];
rz(0.38485011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5220149) q[2];
sx q[2];
rz(-1.3238944) q[2];
sx q[2];
rz(-2.1066693) q[2];
rz(-0.50152913) q[3];
sx q[3];
rz(-0.087567121) q[3];
sx q[3];
rz(-2.1653304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4185249) q[0];
sx q[0];
rz(-1.1534961) q[0];
sx q[0];
rz(-1.204741) q[0];
rz(2.095626) q[1];
sx q[1];
rz(-3.0515262) q[1];
sx q[1];
rz(-3.0012896) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2818272) q[0];
sx q[0];
rz(-0.66376309) q[0];
sx q[0];
rz(1.1361994) q[0];
rz(0.84425462) q[2];
sx q[2];
rz(-1.2469957) q[2];
sx q[2];
rz(2.5716788) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0187776) q[1];
sx q[1];
rz(-0.63624708) q[1];
sx q[1];
rz(-2.754209) q[1];
x q[2];
rz(0.14032538) q[3];
sx q[3];
rz(-1.232339) q[3];
sx q[3];
rz(0.1986157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3673765) q[2];
sx q[2];
rz(-0.99952951) q[2];
sx q[2];
rz(-2.9191169) q[2];
rz(-0.065464822) q[3];
sx q[3];
rz(-1.8583349) q[3];
sx q[3];
rz(-0.84622598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1881994) q[0];
sx q[0];
rz(-3.1173752) q[0];
sx q[0];
rz(-2.5945493) q[0];
rz(-2.8833) q[1];
sx q[1];
rz(-3.1196085) q[1];
sx q[1];
rz(2.8000854) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9753174) q[0];
sx q[0];
rz(-1.5721653) q[0];
sx q[0];
rz(1.5773236) q[0];
rz(1.9893622) q[2];
sx q[2];
rz(-1.8822877) q[2];
sx q[2];
rz(0.77889393) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8402183) q[1];
sx q[1];
rz(-0.84053381) q[1];
sx q[1];
rz(2.6867742) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82587256) q[3];
sx q[3];
rz(-2.4950728) q[3];
sx q[3];
rz(-2.7310128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7121938) q[2];
sx q[2];
rz(-1.2960351) q[2];
sx q[2];
rz(-2.1923501) q[2];
rz(-2.3927169) q[3];
sx q[3];
rz(-1.8604934) q[3];
sx q[3];
rz(3.0794028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6147989) q[0];
sx q[0];
rz(-3.1068046) q[0];
sx q[0];
rz(1.590796) q[0];
rz(1.7855478) q[1];
sx q[1];
rz(-0.0043914774) q[1];
sx q[1];
rz(3.0785676) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7558799) q[0];
sx q[0];
rz(-1.5405419) q[0];
sx q[0];
rz(-0.066019375) q[0];
x q[1];
rz(2.3895415) q[2];
sx q[2];
rz(-1.9785441) q[2];
sx q[2];
rz(-2.4501981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5668727) q[1];
sx q[1];
rz(-1.5441455) q[1];
sx q[1];
rz(1.779056) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0435366) q[3];
sx q[3];
rz(-0.84690988) q[3];
sx q[3];
rz(0.19886097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.65681347) q[2];
sx q[2];
rz(-1.8317089) q[2];
sx q[2];
rz(2.4978034) q[2];
rz(0.8748318) q[3];
sx q[3];
rz(-0.30431408) q[3];
sx q[3];
rz(-0.8005825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0952045) q[0];
sx q[0];
rz(-3.0859741) q[0];
sx q[0];
rz(2.7860506) q[0];
rz(-2.9429759) q[1];
sx q[1];
rz(-0.0067409975) q[1];
sx q[1];
rz(-0.14828646) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3902238) q[0];
sx q[0];
rz(-1.3877739) q[0];
sx q[0];
rz(-1.5682632) q[0];
rz(-1.4951493) q[2];
sx q[2];
rz(-1.9809857) q[2];
sx q[2];
rz(-2.4112521) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1634341) q[1];
sx q[1];
rz(-2.8712007) q[1];
sx q[1];
rz(0.52537523) q[1];
rz(-0.96252302) q[3];
sx q[3];
rz(-1.7136765) q[3];
sx q[3];
rz(2.8382728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6716914) q[2];
sx q[2];
rz(-0.24158676) q[2];
sx q[2];
rz(3.0174603) q[2];
rz(-2.5668674) q[3];
sx q[3];
rz(-0.14437965) q[3];
sx q[3];
rz(2.9706484) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8827051) q[0];
sx q[0];
rz(-0.12508617) q[0];
sx q[0];
rz(2.4001154) q[0];
rz(2.8575836) q[1];
sx q[1];
rz(-0.0037071204) q[1];
sx q[1];
rz(2.8264118) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7595939) q[0];
sx q[0];
rz(-1.6363417) q[0];
sx q[0];
rz(-0.026897341) q[0];
rz(1.1576596) q[2];
sx q[2];
rz(-2.6584932) q[2];
sx q[2];
rz(-0.60483067) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.3005581) q[1];
sx q[1];
rz(-1.4125287) q[1];
sx q[1];
rz(-2.5421418) q[1];
x q[2];
rz(1.4676845) q[3];
sx q[3];
rz(-2.9167843) q[3];
sx q[3];
rz(1.6112427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86889851) q[2];
sx q[2];
rz(-2.0468476) q[2];
sx q[2];
rz(0.72186738) q[2];
rz(2.7591211) q[3];
sx q[3];
rz(-1.1451274) q[3];
sx q[3];
rz(-2.0731879) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.566399) q[0];
sx q[0];
rz(-0.02481758) q[0];
sx q[0];
rz(1.5665293) q[0];
rz(2.9381835) q[1];
sx q[1];
rz(-1.8433488) q[1];
sx q[1];
rz(-2.4967616) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5602556) q[0];
sx q[0];
rz(-1.3632717) q[0];
sx q[0];
rz(0.0085212747) q[0];
rz(0.036083607) q[2];
sx q[2];
rz(-2.5390194) q[2];
sx q[2];
rz(0.35843231) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.825612) q[1];
sx q[1];
rz(-0.54345268) q[1];
sx q[1];
rz(1.4280591) q[1];
rz(-2.0243778) q[3];
sx q[3];
rz(-1.1226113) q[3];
sx q[3];
rz(-2.3635325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5975534) q[2];
sx q[2];
rz(-0.35352239) q[2];
sx q[2];
rz(-0.37975797) q[2];
rz(-1.0455658) q[3];
sx q[3];
rz(-1.9087722) q[3];
sx q[3];
rz(-1.1988962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7757292) q[0];
sx q[0];
rz(-3.1079223) q[0];
sx q[0];
rz(1.3566383) q[0];
rz(0.44048539) q[1];
sx q[1];
rz(-1.0904652) q[1];
sx q[1];
rz(0.7007362) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3627351) q[0];
sx q[0];
rz(-0.88934169) q[0];
sx q[0];
rz(0.89618857) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3323302) q[2];
sx q[2];
rz(-0.70435134) q[2];
sx q[2];
rz(2.4865502) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.060890843) q[1];
sx q[1];
rz(-1.3403088) q[1];
sx q[1];
rz(2.4031765) q[1];
rz(-pi) q[2];
rz(2.7089617) q[3];
sx q[3];
rz(-1.5745161) q[3];
sx q[3];
rz(1.5100117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.35357722) q[2];
sx q[2];
rz(-0.37166301) q[2];
sx q[2];
rz(1.2865944) q[2];
rz(-2.6364117) q[3];
sx q[3];
rz(-0.44613871) q[3];
sx q[3];
rz(-1.9193468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6219567) q[0];
sx q[0];
rz(-3.0918047) q[0];
sx q[0];
rz(-1.5420445) q[0];
rz(-2.385335) q[1];
sx q[1];
rz(-3.1344963) q[1];
sx q[1];
rz(0.33682987) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3339506) q[0];
sx q[0];
rz(-1.5686638) q[0];
sx q[0];
rz(2.8239408) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8518108) q[2];
sx q[2];
rz(-1.769067) q[2];
sx q[2];
rz(-2.2137583) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.67683631) q[1];
sx q[1];
rz(-3.001323) q[1];
sx q[1];
rz(2.2256779) q[1];
rz(-pi) q[2];
rz(-2.1317066) q[3];
sx q[3];
rz(-1.9282544) q[3];
sx q[3];
rz(1.7021029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6280262) q[2];
sx q[2];
rz(-2.1921373) q[2];
sx q[2];
rz(-0.27964082) q[2];
rz(0.63129342) q[3];
sx q[3];
rz(-2.203233) q[3];
sx q[3];
rz(2.5173371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.30017988) q[0];
sx q[0];
rz(-1.5501839) q[0];
sx q[0];
rz(-1.3612904) q[0];
rz(-2.367876) q[1];
sx q[1];
rz(-0.63540375) q[1];
sx q[1];
rz(0.2159963) q[1];
rz(-2.2778289) q[2];
sx q[2];
rz(-2.2939199) q[2];
sx q[2];
rz(-1.6242956) q[2];
rz(-0.37626304) q[3];
sx q[3];
rz(-1.5496764) q[3];
sx q[3];
rz(-1.5839034) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
