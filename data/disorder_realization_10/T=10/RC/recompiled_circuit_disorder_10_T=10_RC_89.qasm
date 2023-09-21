OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(1.1428042) q[0];
sx q[0];
rz(11.354843) q[0];
rz(2.9149574) q[1];
sx q[1];
rz(-1.5645138) q[1];
sx q[1];
rz(-0.29830631) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5814712) q[0];
sx q[0];
rz(-1.6257964) q[0];
sx q[0];
rz(-0.66530692) q[0];
rz(-pi) q[1];
rz(-0.093703336) q[2];
sx q[2];
rz(-1.2393349) q[2];
sx q[2];
rz(-1.9625488) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.28042291) q[1];
sx q[1];
rz(-1.0958332) q[1];
sx q[1];
rz(-0.95649398) q[1];
rz(-pi) q[2];
x q[2];
rz(0.084007752) q[3];
sx q[3];
rz(-1.6735895) q[3];
sx q[3];
rz(1.8889129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6478708) q[2];
sx q[2];
rz(-1.2080668) q[2];
sx q[2];
rz(0.99386627) q[2];
rz(0.99938756) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(-2.5527111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.4988929) q[0];
sx q[0];
rz(-1.9357341) q[0];
sx q[0];
rz(-3.1233741) q[0];
rz(2.3253564) q[1];
sx q[1];
rz(-2.1111592) q[1];
sx q[1];
rz(2.6699064) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52662151) q[0];
sx q[0];
rz(-1.3437628) q[0];
sx q[0];
rz(1.660166) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50498982) q[2];
sx q[2];
rz(-1.8171176) q[2];
sx q[2];
rz(1.1438952) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.43863152) q[1];
sx q[1];
rz(-1.8384117) q[1];
sx q[1];
rz(-2.5665934) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15700335) q[3];
sx q[3];
rz(-2.1727441) q[3];
sx q[3];
rz(-2.4412145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5919684) q[2];
sx q[2];
rz(-1.1529808) q[2];
sx q[2];
rz(-0.96898752) q[2];
rz(0.5747059) q[3];
sx q[3];
rz(-0.55137268) q[3];
sx q[3];
rz(1.0415174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7085003) q[0];
sx q[0];
rz(-1.0555462) q[0];
sx q[0];
rz(-2.95978) q[0];
rz(-1.1026985) q[1];
sx q[1];
rz(-1.6405374) q[1];
sx q[1];
rz(1.6859432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3464976) q[0];
sx q[0];
rz(-0.99612757) q[0];
sx q[0];
rz(-0.57976188) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45708926) q[2];
sx q[2];
rz(-2.3466913) q[2];
sx q[2];
rz(-0.41321102) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0297444) q[1];
sx q[1];
rz(-1.4451471) q[1];
sx q[1];
rz(1.5813584) q[1];
rz(-pi) q[2];
rz(0.62526838) q[3];
sx q[3];
rz(-1.6571952) q[3];
sx q[3];
rz(-2.3141935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2640947) q[2];
sx q[2];
rz(-1.4870746) q[2];
sx q[2];
rz(0.96763119) q[2];
rz(-0.72757059) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(0.23770604) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2847292) q[0];
sx q[0];
rz(-0.52607042) q[0];
sx q[0];
rz(2.5033584) q[0];
rz(2.0137285) q[1];
sx q[1];
rz(-2.3141839) q[1];
sx q[1];
rz(-1.2329873) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5467984) q[0];
sx q[0];
rz(-1.6398755) q[0];
sx q[0];
rz(-0.49537201) q[0];
rz(-2.5523283) q[2];
sx q[2];
rz(-2.0406084) q[2];
sx q[2];
rz(1.1772917) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2265046) q[1];
sx q[1];
rz(-2.9595032) q[1];
sx q[1];
rz(-1.3703129) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2992371) q[3];
sx q[3];
rz(-0.71112594) q[3];
sx q[3];
rz(1.3815051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91810742) q[2];
sx q[2];
rz(-1.2780259) q[2];
sx q[2];
rz(-2.3275862) q[2];
rz(-1.043184) q[3];
sx q[3];
rz(-2.510575) q[3];
sx q[3];
rz(-1.1842747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2994613) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(2.2498851) q[0];
rz(-1.8978329) q[1];
sx q[1];
rz(-1.7638821) q[1];
sx q[1];
rz(0.2125425) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.034886995) q[0];
sx q[0];
rz(-1.5505152) q[0];
sx q[0];
rz(0.01625343) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96617713) q[2];
sx q[2];
rz(-1.539131) q[2];
sx q[2];
rz(2.4633173) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37413874) q[1];
sx q[1];
rz(-0.98857388) q[1];
sx q[1];
rz(0.9064845) q[1];
rz(-2.2977923) q[3];
sx q[3];
rz(-1.926933) q[3];
sx q[3];
rz(2.2775377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.033096878) q[2];
sx q[2];
rz(-0.97110811) q[2];
sx q[2];
rz(0.5212211) q[2];
rz(1.7565578) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(-1.8732171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2824771) q[0];
sx q[0];
rz(-0.92882597) q[0];
sx q[0];
rz(-2.916472) q[0];
rz(-1.7865932) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(0.37757847) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7901944) q[0];
sx q[0];
rz(-1.5488008) q[0];
sx q[0];
rz(-1.8368594) q[0];
rz(-0.2874561) q[2];
sx q[2];
rz(-1.6171347) q[2];
sx q[2];
rz(-3.1099144) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0203591) q[1];
sx q[1];
rz(-1.6755783) q[1];
sx q[1];
rz(1.1377513) q[1];
x q[2];
rz(-0.79490957) q[3];
sx q[3];
rz(-1.5494293) q[3];
sx q[3];
rz(1.6705318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.98465115) q[2];
sx q[2];
rz(-2.3539383) q[2];
sx q[2];
rz(2.6605576) q[2];
rz(-2.7379819) q[3];
sx q[3];
rz(-2.1026881) q[3];
sx q[3];
rz(2.8267982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0525381) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(2.9470434) q[0];
rz(0.21952195) q[1];
sx q[1];
rz(-1.4621282) q[1];
sx q[1];
rz(0.25442466) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62691488) q[0];
sx q[0];
rz(-1.2509545) q[0];
sx q[0];
rz(0.11492782) q[0];
rz(2.9872586) q[2];
sx q[2];
rz(-0.96417226) q[2];
sx q[2];
rz(2.2301205) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.82032694) q[1];
sx q[1];
rz(-1.1202381) q[1];
sx q[1];
rz(1.2688368) q[1];
x q[2];
rz(-2.3354704) q[3];
sx q[3];
rz(-0.61180173) q[3];
sx q[3];
rz(2.7968189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.97757942) q[2];
sx q[2];
rz(-1.7501202) q[2];
sx q[2];
rz(-1.3158201) q[2];
rz(0.87604648) q[3];
sx q[3];
rz(-0.13893572) q[3];
sx q[3];
rz(-1.0036489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2106237) q[0];
sx q[0];
rz(-0.3759149) q[0];
sx q[0];
rz(1.6865431) q[0];
rz(2.3176106) q[1];
sx q[1];
rz(-2.1897557) q[1];
sx q[1];
rz(-1.5751858) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63021916) q[0];
sx q[0];
rz(-0.89238088) q[0];
sx q[0];
rz(1.9766962) q[0];
rz(-2.1112061) q[2];
sx q[2];
rz(-1.4166797) q[2];
sx q[2];
rz(2.266778) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34711429) q[1];
sx q[1];
rz(-1.9273259) q[1];
sx q[1];
rz(-1.1785248) q[1];
rz(-pi) q[2];
rz(-0.67993645) q[3];
sx q[3];
rz(-0.61938647) q[3];
sx q[3];
rz(-3.1114651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2660797) q[2];
sx q[2];
rz(-1.2337039) q[2];
sx q[2];
rz(2.8640462) q[2];
rz(-1.6905486) q[3];
sx q[3];
rz(-2.6896559) q[3];
sx q[3];
rz(-1.014876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16185109) q[0];
sx q[0];
rz(-2.6888872) q[0];
sx q[0];
rz(-1.45654) q[0];
rz(2.5121571) q[1];
sx q[1];
rz(-1.9742191) q[1];
sx q[1];
rz(1.1368407) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0757383) q[0];
sx q[0];
rz(-0.59959164) q[0];
sx q[0];
rz(-1.0435186) q[0];
rz(-2.7669737) q[2];
sx q[2];
rz(-1.9667452) q[2];
sx q[2];
rz(-1.7916726) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7965664) q[1];
sx q[1];
rz(-1.476164) q[1];
sx q[1];
rz(1.8887397) q[1];
rz(1.198248) q[3];
sx q[3];
rz(-2.2714943) q[3];
sx q[3];
rz(-2.8230599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64951605) q[2];
sx q[2];
rz(-1.7121544) q[2];
sx q[2];
rz(-2.1006404) q[2];
rz(3.1395636) q[3];
sx q[3];
rz(-2.7414331) q[3];
sx q[3];
rz(2.8295529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8390389) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(2.1955406) q[0];
rz(0.91167766) q[1];
sx q[1];
rz(-1.2152351) q[1];
sx q[1];
rz(-0.61202234) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34110585) q[0];
sx q[0];
rz(-2.2414811) q[0];
sx q[0];
rz(0.97408803) q[0];
x q[1];
rz(-2.5435796) q[2];
sx q[2];
rz(-0.41053718) q[2];
sx q[2];
rz(1.7288127) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.40572383) q[1];
sx q[1];
rz(-1.6940261) q[1];
sx q[1];
rz(0.65995364) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53491433) q[3];
sx q[3];
rz(-1.8992918) q[3];
sx q[3];
rz(1.2319777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0570021) q[2];
sx q[2];
rz(-2.4995063) q[2];
sx q[2];
rz(2.4882312) q[2];
rz(-0.35081321) q[3];
sx q[3];
rz(-1.6143129) q[3];
sx q[3];
rz(-2.4408834) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54031298) q[0];
sx q[0];
rz(-1.5355587) q[0];
sx q[0];
rz(-3.0328947) q[0];
rz(0.75469771) q[1];
sx q[1];
rz(-1.8221868) q[1];
sx q[1];
rz(1.6356161) q[1];
rz(-2.4621261) q[2];
sx q[2];
rz(-1.1107399) q[2];
sx q[2];
rz(-2.6723292) q[2];
rz(-3.0572206) q[3];
sx q[3];
rz(-1.2481239) q[3];
sx q[3];
rz(0.39831755) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];