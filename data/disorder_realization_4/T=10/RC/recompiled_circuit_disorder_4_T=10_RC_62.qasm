OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274347) q[0];
sx q[0];
rz(-0.43570575) q[0];
sx q[0];
rz(0.92619196) q[0];
rz(1.9594833) q[1];
sx q[1];
rz(-0.73298454) q[1];
sx q[1];
rz(-2.7690673) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9076219) q[0];
sx q[0];
rz(-1.3682433) q[0];
sx q[0];
rz(-2.5552208) q[0];
x q[1];
rz(1.4650605) q[2];
sx q[2];
rz(-0.9423965) q[2];
sx q[2];
rz(2.489593) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.26973154) q[1];
sx q[1];
rz(-1.6850123) q[1];
sx q[1];
rz(-1.7608587) q[1];
rz(3.0285809) q[3];
sx q[3];
rz(-1.3556619) q[3];
sx q[3];
rz(0.61258951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7131876) q[2];
sx q[2];
rz(-1.5822072) q[2];
sx q[2];
rz(-0.92450809) q[2];
rz(-1.6690147) q[3];
sx q[3];
rz(-1.2481097) q[3];
sx q[3];
rz(-1.4991466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(2.9532042) q[0];
sx q[0];
rz(-0.062347118) q[0];
sx q[0];
rz(-1.6054608) q[0];
rz(-0.19451441) q[1];
sx q[1];
rz(-1.3214) q[1];
sx q[1];
rz(3.0867192) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7354436) q[0];
sx q[0];
rz(-2.9950954) q[0];
sx q[0];
rz(-2.2551401) q[0];
rz(-pi) q[1];
rz(-0.70116455) q[2];
sx q[2];
rz(-0.75116457) q[2];
sx q[2];
rz(-3.1322111) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4108737) q[1];
sx q[1];
rz(-1.5667856) q[1];
sx q[1];
rz(1.2282759) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7826471) q[3];
sx q[3];
rz(-2.2329674) q[3];
sx q[3];
rz(0.3609095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7028971) q[2];
sx q[2];
rz(-0.56151152) q[2];
sx q[2];
rz(-1.4820209) q[2];
rz(2.1510018) q[3];
sx q[3];
rz(-1.7329268) q[3];
sx q[3];
rz(-1.3668485) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3593339) q[0];
sx q[0];
rz(-2.6222836) q[0];
sx q[0];
rz(-2.7666132) q[0];
rz(0.24770501) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(-1.3365655) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61993116) q[0];
sx q[0];
rz(-1.7921899) q[0];
sx q[0];
rz(3.0822166) q[0];
rz(-pi) q[1];
rz(0.61632421) q[2];
sx q[2];
rz(-2.4547572) q[2];
sx q[2];
rz(0.64885215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.49866906) q[1];
sx q[1];
rz(-1.8734697) q[1];
sx q[1];
rz(-2.32294) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1683447) q[3];
sx q[3];
rz(-2.2964722) q[3];
sx q[3];
rz(-2.890051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1027362) q[2];
sx q[2];
rz(-2.3114624) q[2];
sx q[2];
rz(-2.1170199) q[2];
rz(1.2767977) q[3];
sx q[3];
rz(-1.5921311) q[3];
sx q[3];
rz(-1.6259441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643726) q[0];
sx q[0];
rz(-1.6765046) q[0];
sx q[0];
rz(-3.1080416) q[0];
rz(1.230348) q[1];
sx q[1];
rz(-0.76554811) q[1];
sx q[1];
rz(-1.4039325) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0516292) q[0];
sx q[0];
rz(-1.3768059) q[0];
sx q[0];
rz(0.35332638) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83909859) q[2];
sx q[2];
rz(-2.1144146) q[2];
sx q[2];
rz(0.77857882) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0714604) q[1];
sx q[1];
rz(-1.4590108) q[1];
sx q[1];
rz(1.2791355) q[1];
x q[2];
rz(2.2622044) q[3];
sx q[3];
rz(-1.6367568) q[3];
sx q[3];
rz(-1.1018167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6491062) q[2];
sx q[2];
rz(-2.2968569) q[2];
sx q[2];
rz(0.21326324) q[2];
rz(-0.02031859) q[3];
sx q[3];
rz(-1.820194) q[3];
sx q[3];
rz(0.4030574) q[3];
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
rz(-1.7810818) q[0];
sx q[0];
rz(-1.1163982) q[0];
sx q[0];
rz(2.1257341) q[0];
rz(1.6332743) q[1];
sx q[1];
rz(-0.95817482) q[1];
sx q[1];
rz(3.1075081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.236892) q[0];
sx q[0];
rz(-1.1780945) q[0];
sx q[0];
rz(-1.3971726) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57406117) q[2];
sx q[2];
rz(-1.7957557) q[2];
sx q[2];
rz(-3.0846734) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2195697) q[1];
sx q[1];
rz(-1.2864188) q[1];
sx q[1];
rz(-0.0071830458) q[1];
rz(0.91655101) q[3];
sx q[3];
rz(-0.74186462) q[3];
sx q[3];
rz(1.3548917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.74026996) q[2];
sx q[2];
rz(-2.5677742) q[2];
sx q[2];
rz(1.1711228) q[2];
rz(-1.8088388) q[3];
sx q[3];
rz(-1.4274024) q[3];
sx q[3];
rz(-0.12935054) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68833441) q[0];
sx q[0];
rz(-1.9535221) q[0];
sx q[0];
rz(-2.7057498) q[0];
rz(-2.0360937) q[1];
sx q[1];
rz(-1.8445797) q[1];
sx q[1];
rz(1.2058535) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0379369) q[0];
sx q[0];
rz(-2.3559542) q[0];
sx q[0];
rz(1.2334137) q[0];
rz(-pi) q[1];
rz(0.016024307) q[2];
sx q[2];
rz(-1.4756225) q[2];
sx q[2];
rz(-1.3118088) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1356493) q[1];
sx q[1];
rz(-1.5444396) q[1];
sx q[1];
rz(-2.6308358) q[1];
rz(-pi) q[2];
rz(1.1711575) q[3];
sx q[3];
rz(-1.6013147) q[3];
sx q[3];
rz(0.37351028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2445406) q[2];
sx q[2];
rz(-2.5138833) q[2];
sx q[2];
rz(-0.051606027) q[2];
rz(-0.40766019) q[3];
sx q[3];
rz(-1.1838653) q[3];
sx q[3];
rz(0.44529644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.62149858) q[0];
sx q[0];
rz(-0.61722732) q[0];
sx q[0];
rz(0.67333418) q[0];
rz(-2.3576221) q[1];
sx q[1];
rz(-1.7242804) q[1];
sx q[1];
rz(0.53692445) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2818031) q[0];
sx q[0];
rz(-1.4872695) q[0];
sx q[0];
rz(-3.0755755) q[0];
x q[1];
rz(-1.827048) q[2];
sx q[2];
rz(-1.7125687) q[2];
sx q[2];
rz(-1.2830551) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6785477) q[1];
sx q[1];
rz(-1.6650513) q[1];
sx q[1];
rz(-2.3386392) q[1];
rz(-pi) q[2];
rz(0.22646871) q[3];
sx q[3];
rz(-2.5385751) q[3];
sx q[3];
rz(-0.23803593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15988222) q[2];
sx q[2];
rz(-1.8843702) q[2];
sx q[2];
rz(2.9505777) q[2];
rz(2.8602709) q[3];
sx q[3];
rz(-1.9502935) q[3];
sx q[3];
rz(-0.34240001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1087588) q[0];
sx q[0];
rz(-0.37029752) q[0];
sx q[0];
rz(-2.7539745) q[0];
rz(3.0265813) q[1];
sx q[1];
rz(-1.4287881) q[1];
sx q[1];
rz(0.33755916) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.176929) q[0];
sx q[0];
rz(-1.3725855) q[0];
sx q[0];
rz(0.49074178) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6667716) q[2];
sx q[2];
rz(-2.5002067) q[2];
sx q[2];
rz(2.1542187) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.90606373) q[1];
sx q[1];
rz(-1.4396832) q[1];
sx q[1];
rz(-2.1436585) q[1];
rz(-pi) q[2];
rz(0.75546219) q[3];
sx q[3];
rz(-1.6172234) q[3];
sx q[3];
rz(-1.69343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.18743029) q[2];
sx q[2];
rz(-0.53415161) q[2];
sx q[2];
rz(1.2672651) q[2];
rz(0.77504843) q[3];
sx q[3];
rz(-1.5379484) q[3];
sx q[3];
rz(2.3601941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18701126) q[0];
sx q[0];
rz(-1.0382074) q[0];
sx q[0];
rz(2.9206081) q[0];
rz(-2.2194608) q[1];
sx q[1];
rz(-1.2698959) q[1];
sx q[1];
rz(-1.0029213) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98175883) q[0];
sx q[0];
rz(-0.28875414) q[0];
sx q[0];
rz(-0.50555484) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1557751) q[2];
sx q[2];
rz(-1.1406116) q[2];
sx q[2];
rz(-2.1625724) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6332615) q[1];
sx q[1];
rz(-1.08053) q[1];
sx q[1];
rz(-1.0432748) q[1];
rz(-pi) q[2];
rz(-0.88463155) q[3];
sx q[3];
rz(-1.7836708) q[3];
sx q[3];
rz(-0.5205982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6691436) q[2];
sx q[2];
rz(-2.3995235) q[2];
sx q[2];
rz(-0.15850244) q[2];
rz(-0.45378271) q[3];
sx q[3];
rz(-2.3588389) q[3];
sx q[3];
rz(-0.84428549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
rz(2.6145988) q[0];
sx q[0];
rz(-2.232593) q[0];
sx q[0];
rz(0.19113834) q[0];
rz(2.846431) q[1];
sx q[1];
rz(-2.252153) q[1];
sx q[1];
rz(-2.2492762) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79217171) q[0];
sx q[0];
rz(-0.76602174) q[0];
sx q[0];
rz(1.8269405) q[0];
x q[1];
rz(1.4453663) q[2];
sx q[2];
rz(-2.4969366) q[2];
sx q[2];
rz(0.21412011) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.416623) q[1];
sx q[1];
rz(-1.6367216) q[1];
sx q[1];
rz(2.9030187) q[1];
rz(-pi) q[2];
rz(0.43182208) q[3];
sx q[3];
rz(-1.338827) q[3];
sx q[3];
rz(-0.98917978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7717379) q[2];
sx q[2];
rz(-2.3042046) q[2];
sx q[2];
rz(0.85956335) q[2];
rz(-1.9178948) q[3];
sx q[3];
rz(-0.92646354) q[3];
sx q[3];
rz(2.3971476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0201465) q[0];
sx q[0];
rz(-0.93562026) q[0];
sx q[0];
rz(1.390441) q[0];
rz(1.4795115) q[1];
sx q[1];
rz(-0.48304396) q[1];
sx q[1];
rz(1.2190291) q[1];
rz(0.027364846) q[2];
sx q[2];
rz(-1.8623427) q[2];
sx q[2];
rz(1.3775415) q[2];
rz(1.6395232) q[3];
sx q[3];
rz(-1.6106265) q[3];
sx q[3];
rz(2.0911218) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];