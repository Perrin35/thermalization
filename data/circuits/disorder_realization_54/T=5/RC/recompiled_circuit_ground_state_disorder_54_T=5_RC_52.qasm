OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7689826) q[0];
sx q[0];
rz(3.1867653) q[0];
sx q[0];
rz(10.09633) q[0];
rz(2.1454732) q[1];
sx q[1];
rz(5.4944333) q[1];
sx q[1];
rz(6.5715437) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2644293) q[0];
sx q[0];
rz(-0.20119871) q[0];
sx q[0];
rz(1.3876794) q[0];
rz(-pi) q[1];
rz(1.8419349) q[2];
sx q[2];
rz(-2.2170005) q[2];
sx q[2];
rz(0.14474511) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3934196) q[1];
sx q[1];
rz(-0.44965023) q[1];
sx q[1];
rz(-0.086918615) q[1];
x q[2];
rz(0.45436556) q[3];
sx q[3];
rz(-1.7011257) q[3];
sx q[3];
rz(-0.59727515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.70197376) q[2];
sx q[2];
rz(-0.34586033) q[2];
sx q[2];
rz(2.4227179) q[2];
rz(-1.4465205) q[3];
sx q[3];
rz(-1.6547763) q[3];
sx q[3];
rz(-1.0369302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.059939) q[0];
sx q[0];
rz(-0.43123284) q[0];
sx q[0];
rz(2.8531139) q[0];
rz(2.5892995) q[1];
sx q[1];
rz(-1.0918795) q[1];
sx q[1];
rz(1.1757895) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8771851) q[0];
sx q[0];
rz(-2.7938936) q[0];
sx q[0];
rz(-2.023979) q[0];
x q[1];
rz(-1.7442877) q[2];
sx q[2];
rz(-1.6855006) q[2];
sx q[2];
rz(-0.8300654) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.77681357) q[1];
sx q[1];
rz(-2.5506667) q[1];
sx q[1];
rz(1.7950115) q[1];
rz(-pi) q[2];
rz(1.5351686) q[3];
sx q[3];
rz(-3.1094915) q[3];
sx q[3];
rz(2.7184582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6931849) q[2];
sx q[2];
rz(-1.2625445) q[2];
sx q[2];
rz(0.022424879) q[2];
rz(1.4261931) q[3];
sx q[3];
rz(-1.8636999) q[3];
sx q[3];
rz(2.2414331) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65396032) q[0];
sx q[0];
rz(-1.2314236) q[0];
sx q[0];
rz(-2.3768429) q[0];
rz(-1.7817616) q[1];
sx q[1];
rz(-1.9197074) q[1];
sx q[1];
rz(-1.2545895) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6343289) q[0];
sx q[0];
rz(-1.467469) q[0];
sx q[0];
rz(1.4402585) q[0];
rz(-pi) q[1];
rz(-2.8171982) q[2];
sx q[2];
rz(-1.430871) q[2];
sx q[2];
rz(3.0382699) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0618825) q[1];
sx q[1];
rz(-1.2822423) q[1];
sx q[1];
rz(0.48019073) q[1];
rz(-pi) q[2];
rz(-1.4419008) q[3];
sx q[3];
rz(-1.2308972) q[3];
sx q[3];
rz(1.631402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7356871) q[2];
sx q[2];
rz(-0.76852208) q[2];
sx q[2];
rz(-0.020966919) q[2];
rz(1.5603125) q[3];
sx q[3];
rz(-2.0798648) q[3];
sx q[3];
rz(2.235152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2727994) q[0];
sx q[0];
rz(-0.11196207) q[0];
sx q[0];
rz(1.3695166) q[0];
rz(0.70961332) q[1];
sx q[1];
rz(-1.2779002) q[1];
sx q[1];
rz(2.4443464) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.227237) q[0];
sx q[0];
rz(-0.52292693) q[0];
sx q[0];
rz(0.65547734) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68402779) q[2];
sx q[2];
rz(-2.1993756) q[2];
sx q[2];
rz(-2.7492439) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5124951) q[1];
sx q[1];
rz(-1.2689018) q[1];
sx q[1];
rz(-0.90934335) q[1];
rz(-0.94631291) q[3];
sx q[3];
rz(-0.28879228) q[3];
sx q[3];
rz(-0.88209796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9559481) q[2];
sx q[2];
rz(-2.1944025) q[2];
sx q[2];
rz(0.67898018) q[2];
rz(-1.125157) q[3];
sx q[3];
rz(-1.1449287) q[3];
sx q[3];
rz(0.2230491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0515902) q[0];
sx q[0];
rz(-1.2035878) q[0];
sx q[0];
rz(-3.0644655) q[0];
rz(-1.1513101) q[1];
sx q[1];
rz(-1.5733893) q[1];
sx q[1];
rz(0.77879771) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1098233) q[0];
sx q[0];
rz(-1.6187877) q[0];
sx q[0];
rz(0.035295156) q[0];
rz(1.2856712) q[2];
sx q[2];
rz(-1.4094009) q[2];
sx q[2];
rz(1.6708167) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.59757876) q[1];
sx q[1];
rz(-2.72157) q[1];
sx q[1];
rz(2.2060664) q[1];
rz(-pi) q[2];
rz(-0.058249931) q[3];
sx q[3];
rz(-0.16213972) q[3];
sx q[3];
rz(-1.8193965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6246346) q[2];
sx q[2];
rz(-1.4411074) q[2];
sx q[2];
rz(-1.9514294) q[2];
rz(0.55365753) q[3];
sx q[3];
rz(-2.5167969) q[3];
sx q[3];
rz(-2.4042118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6051642) q[0];
sx q[0];
rz(-1.9909415) q[0];
sx q[0];
rz(-2.8705257) q[0];
rz(-1.0003264) q[1];
sx q[1];
rz(-1.9258291) q[1];
sx q[1];
rz(-2.2727374) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1329089) q[0];
sx q[0];
rz(-2.0374679) q[0];
sx q[0];
rz(-2.4746289) q[0];
rz(2.4060288) q[2];
sx q[2];
rz(-1.6297139) q[2];
sx q[2];
rz(-1.9933188) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0583494) q[1];
sx q[1];
rz(-1.740137) q[1];
sx q[1];
rz(1.9237299) q[1];
rz(-pi) q[2];
rz(-2.3575918) q[3];
sx q[3];
rz(-1.2262359) q[3];
sx q[3];
rz(1.6387516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.102999) q[2];
sx q[2];
rz(-2.8391892) q[2];
sx q[2];
rz(-2.0043066) q[2];
rz(0.31050995) q[3];
sx q[3];
rz(-2.2313084) q[3];
sx q[3];
rz(2.212132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1167269) q[0];
sx q[0];
rz(-1.4288582) q[0];
sx q[0];
rz(-2.0615935) q[0];
rz(-2.179821) q[1];
sx q[1];
rz(-0.56517833) q[1];
sx q[1];
rz(-2.9186509) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38755435) q[0];
sx q[0];
rz(-1.5359274) q[0];
sx q[0];
rz(-0.0074444093) q[0];
x q[1];
rz(2.8330363) q[2];
sx q[2];
rz(-2.1962104) q[2];
sx q[2];
rz(2.5297414) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2126951) q[1];
sx q[1];
rz(-1.2948117) q[1];
sx q[1];
rz(2.3174965) q[1];
rz(0.014752905) q[3];
sx q[3];
rz(-0.97501576) q[3];
sx q[3];
rz(2.9500913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.96233931) q[2];
sx q[2];
rz(-0.27955678) q[2];
sx q[2];
rz(-1.9635828) q[2];
rz(2.8152605) q[3];
sx q[3];
rz(-0.81063619) q[3];
sx q[3];
rz(-2.0012205) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2278263) q[0];
sx q[0];
rz(-0.95781177) q[0];
sx q[0];
rz(0.80818278) q[0];
rz(0.35762865) q[1];
sx q[1];
rz(-1.459815) q[1];
sx q[1];
rz(-1.0825895) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4935166) q[0];
sx q[0];
rz(-2.2449623) q[0];
sx q[0];
rz(-0.21866673) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0183187) q[2];
sx q[2];
rz(-1.7874771) q[2];
sx q[2];
rz(-2.0683757) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.030176) q[1];
sx q[1];
rz(-1.7716265) q[1];
sx q[1];
rz(-2.3371731) q[1];
rz(-2.0478422) q[3];
sx q[3];
rz(-1.1594605) q[3];
sx q[3];
rz(-0.75440948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0766582) q[2];
sx q[2];
rz(-0.4883464) q[2];
sx q[2];
rz(-1.2154382) q[2];
rz(-0.91228929) q[3];
sx q[3];
rz(-0.89022294) q[3];
sx q[3];
rz(0.54273763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2927581) q[0];
sx q[0];
rz(-0.64105761) q[0];
sx q[0];
rz(2.7178398) q[0];
rz(-1.1774225) q[1];
sx q[1];
rz(-2.2260428) q[1];
sx q[1];
rz(-0.37568572) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59058023) q[0];
sx q[0];
rz(-2.9504021) q[0];
sx q[0];
rz(1.1893597) q[0];
rz(-pi) q[1];
x q[1];
rz(3.032349) q[2];
sx q[2];
rz(-2.901361) q[2];
sx q[2];
rz(1.9640528) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0547395) q[1];
sx q[1];
rz(-2.6356831) q[1];
sx q[1];
rz(1.4356126) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71640941) q[3];
sx q[3];
rz(-1.3103974) q[3];
sx q[3];
rz(1.474787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11289135) q[2];
sx q[2];
rz(-1.4914923) q[2];
sx q[2];
rz(2.4493307) q[2];
rz(-0.68474692) q[3];
sx q[3];
rz(-2.1919577) q[3];
sx q[3];
rz(-0.14028604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.58458644) q[0];
sx q[0];
rz(-0.96730119) q[0];
sx q[0];
rz(-1.6242356) q[0];
rz(1.5380305) q[1];
sx q[1];
rz(-1.3471194) q[1];
sx q[1];
rz(1.2204407) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6950078) q[0];
sx q[0];
rz(-1.5112226) q[0];
sx q[0];
rz(2.1936962) q[0];
rz(-0.50165711) q[2];
sx q[2];
rz(-1.3634472) q[2];
sx q[2];
rz(-0.32688552) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5511849) q[1];
sx q[1];
rz(-0.77279323) q[1];
sx q[1];
rz(2.7332952) q[1];
rz(-pi) q[2];
rz(2.4046005) q[3];
sx q[3];
rz(-1.3496282) q[3];
sx q[3];
rz(-1.8870518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5662235) q[2];
sx q[2];
rz(-1.0074002) q[2];
sx q[2];
rz(-2.2929906) q[2];
rz(-0.97950116) q[3];
sx q[3];
rz(-1.5749911) q[3];
sx q[3];
rz(0.037467329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.033087) q[0];
sx q[0];
rz(-1.6642234) q[0];
sx q[0];
rz(-2.4432175) q[0];
rz(-2.5820844) q[1];
sx q[1];
rz(-2.5026176) q[1];
sx q[1];
rz(2.7284596) q[1];
rz(-2.8200061) q[2];
sx q[2];
rz(-1.3272616) q[2];
sx q[2];
rz(2.7250901) q[2];
rz(-1.63089) q[3];
sx q[3];
rz(-0.99219764) q[3];
sx q[3];
rz(0.037430684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
