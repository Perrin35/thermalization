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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059271185) q[0];
sx q[0];
rz(-2.4743609) q[0];
sx q[0];
rz(3.0526403) q[0];
x q[1];
rz(1.9036129) q[2];
sx q[2];
rz(-1.6593854) q[2];
sx q[2];
rz(0.36117902) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6023941) q[1];
sx q[1];
rz(-2.1089923) q[1];
sx q[1];
rz(-2.5799275) q[1];
x q[2];
rz(-3.0575849) q[3];
sx q[3];
rz(-1.4680032) q[3];
sx q[3];
rz(-1.8889129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4937218) q[2];
sx q[2];
rz(-1.9335258) q[2];
sx q[2];
rz(-0.99386627) q[2];
rz(2.1422051) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(2.5527111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64269972) q[0];
sx q[0];
rz(-1.2058586) q[0];
sx q[0];
rz(-0.018218536) q[0];
rz(-0.81623626) q[1];
sx q[1];
rz(-2.1111592) q[1];
sx q[1];
rz(-0.47168628) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2361006) q[0];
sx q[0];
rz(-0.24370757) q[0];
sx q[0];
rz(-2.7729176) q[0];
rz(-pi) q[1];
rz(2.6366028) q[2];
sx q[2];
rz(-1.8171176) q[2];
sx q[2];
rz(-1.9976975) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7029611) q[1];
sx q[1];
rz(-1.8384117) q[1];
sx q[1];
rz(-0.57499927) q[1];
x q[2];
rz(1.3470115) q[3];
sx q[3];
rz(-2.5219678) q[3];
sx q[3];
rz(0.97298813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5919684) q[2];
sx q[2];
rz(-1.1529808) q[2];
sx q[2];
rz(2.1726051) q[2];
rz(-2.5668868) q[3];
sx q[3];
rz(-2.59022) q[3];
sx q[3];
rz(-1.0415174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4330924) q[0];
sx q[0];
rz(-2.0860465) q[0];
sx q[0];
rz(2.95978) q[0];
rz(-1.1026985) q[1];
sx q[1];
rz(-1.6405374) q[1];
sx q[1];
rz(1.6859432) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91711125) q[0];
sx q[0];
rz(-2.3492976) q[0];
sx q[0];
rz(-0.86865058) q[0];
x q[1];
rz(1.9934898) q[2];
sx q[2];
rz(-2.2659677) q[2];
sx q[2];
rz(-0.19902755) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.54237565) q[1];
sx q[1];
rz(-1.5603175) q[1];
sx q[1];
rz(0.12565617) q[1];
rz(1.6772179) q[3];
sx q[3];
rz(-0.94821804) q[3];
sx q[3];
rz(2.4604083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87749798) q[2];
sx q[2];
rz(-1.4870746) q[2];
sx q[2];
rz(2.1739615) q[2];
rz(-0.72757059) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(-2.9038866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.85686344) q[0];
sx q[0];
rz(-2.6155222) q[0];
sx q[0];
rz(0.63823429) q[0];
rz(-2.0137285) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(1.9086054) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9905332) q[0];
sx q[0];
rz(-0.49976832) q[0];
sx q[0];
rz(-0.14453669) q[0];
x q[1];
rz(0.74027503) q[2];
sx q[2];
rz(-0.7358272) q[2];
sx q[2];
rz(-2.1528113) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4302664) q[1];
sx q[1];
rz(-1.749199) q[1];
sx q[1];
rz(3.1049411) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2992371) q[3];
sx q[3];
rz(-2.4304667) q[3];
sx q[3];
rz(-1.3815051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2234852) q[2];
sx q[2];
rz(-1.8635668) q[2];
sx q[2];
rz(-0.81400648) q[2];
rz(1.043184) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(-1.1842747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.84213132) q[0];
sx q[0];
rz(-1.786754) q[0];
sx q[0];
rz(0.89170757) q[0];
rz(-1.2437598) q[1];
sx q[1];
rz(-1.3777106) q[1];
sx q[1];
rz(-2.9290501) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71056238) q[0];
sx q[0];
rz(-3.115603) q[0];
sx q[0];
rz(-2.2463069) q[0];
rz(1.5151305) q[2];
sx q[2];
rz(-2.5362483) q[2];
sx q[2];
rz(-0.84673131) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.58442851) q[1];
sx q[1];
rz(-2.2884528) q[1];
sx q[1];
rz(-0.75259705) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6796474) q[3];
sx q[3];
rz(-0.89832234) q[3];
sx q[3];
rz(0.40601054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.033096878) q[2];
sx q[2];
rz(-2.1704845) q[2];
sx q[2];
rz(2.6203716) q[2];
rz(1.7565578) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(1.2683755) q[3];
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
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8591156) q[0];
sx q[0];
rz(-0.92882597) q[0];
sx q[0];
rz(2.916472) q[0];
rz(1.7865932) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(2.7640142) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7901944) q[0];
sx q[0];
rz(-1.5927918) q[0];
sx q[0];
rz(-1.3047332) q[0];
rz(0.2874561) q[2];
sx q[2];
rz(-1.6171347) q[2];
sx q[2];
rz(3.1099144) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4695417) q[1];
sx q[1];
rz(-0.44476032) q[1];
sx q[1];
rz(-1.3252392) q[1];
rz(-2.3466831) q[3];
sx q[3];
rz(-1.5921633) q[3];
sx q[3];
rz(-1.4710609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98465115) q[2];
sx q[2];
rz(-2.3539383) q[2];
sx q[2];
rz(-0.48103508) q[2];
rz(0.40361079) q[3];
sx q[3];
rz(-1.0389046) q[3];
sx q[3];
rz(-2.8267982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0890546) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(-0.19454923) q[0];
rz(0.21952195) q[1];
sx q[1];
rz(-1.4621282) q[1];
sx q[1];
rz(0.25442466) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2339904) q[0];
sx q[0];
rz(-1.6798717) q[0];
sx q[0];
rz(1.8926189) q[0];
rz(-1.3527855) q[2];
sx q[2];
rz(-2.5180452) q[2];
sx q[2];
rz(0.64507285) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3212657) q[1];
sx q[1];
rz(-2.0213545) q[1];
sx q[1];
rz(-1.8727559) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3354704) q[3];
sx q[3];
rz(-0.61180173) q[3];
sx q[3];
rz(2.7968189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.97757942) q[2];
sx q[2];
rz(-1.3914725) q[2];
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
rz(pi/2) q[1];
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
rz(-1.9309689) q[0];
sx q[0];
rz(-2.7656778) q[0];
sx q[0];
rz(1.6865431) q[0];
rz(-0.82398206) q[1];
sx q[1];
rz(-0.95183698) q[1];
sx q[1];
rz(1.5751858) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5113735) q[0];
sx q[0];
rz(-2.2492118) q[0];
sx q[0];
rz(1.9766962) q[0];
x q[1];
rz(1.8640395) q[2];
sx q[2];
rz(-2.5817421) q[2];
sx q[2];
rz(0.44550371) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3670866) q[1];
sx q[1];
rz(-1.204406) q[1];
sx q[1];
rz(0.3831425) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67993645) q[3];
sx q[3];
rz(-2.5222062) q[3];
sx q[3];
rz(-0.030127545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87551293) q[2];
sx q[2];
rz(-1.9078887) q[2];
sx q[2];
rz(-2.8640462) q[2];
rz(-1.4510441) q[3];
sx q[3];
rz(-0.45193672) q[3];
sx q[3];
rz(2.1267166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9797416) q[0];
sx q[0];
rz(-0.45270544) q[0];
sx q[0];
rz(-1.6850527) q[0];
rz(-0.62943554) q[1];
sx q[1];
rz(-1.1673735) q[1];
sx q[1];
rz(2.004752) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6900078) q[0];
sx q[0];
rz(-1.0613872) q[0];
sx q[0];
rz(2.8103229) q[0];
rz(2.7669737) q[2];
sx q[2];
rz(-1.1748474) q[2];
sx q[2];
rz(-1.7916726) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.34502621) q[1];
sx q[1];
rz(-1.476164) q[1];
sx q[1];
rz(1.252853) q[1];
x q[2];
rz(-1.9433446) q[3];
sx q[3];
rz(-0.87009831) q[3];
sx q[3];
rz(2.8230599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.64951605) q[2];
sx q[2];
rz(-1.4294383) q[2];
sx q[2];
rz(2.1006404) q[2];
rz(-3.1395636) q[3];
sx q[3];
rz(-2.7414331) q[3];
sx q[3];
rz(0.31203976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3025538) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(-0.94605207) q[0];
rz(-0.91167766) q[1];
sx q[1];
rz(-1.9263575) q[1];
sx q[1];
rz(-0.61202234) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6291954) q[0];
sx q[0];
rz(-1.1149659) q[0];
sx q[0];
rz(2.3770611) q[0];
rz(-1.811118) q[2];
sx q[2];
rz(-1.2346621) q[2];
sx q[2];
rz(-2.0517595) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7358688) q[1];
sx q[1];
rz(-1.6940261) q[1];
sx q[1];
rz(2.481639) q[1];
rz(-pi) q[2];
rz(-2.5522334) q[3];
sx q[3];
rz(-0.61925626) q[3];
sx q[3];
rz(-2.9818231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0570021) q[2];
sx q[2];
rz(-0.64208639) q[2];
sx q[2];
rz(2.4882312) q[2];
rz(2.7907794) q[3];
sx q[3];
rz(-1.6143129) q[3];
sx q[3];
rz(0.70070926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6012797) q[0];
sx q[0];
rz(-1.606034) q[0];
sx q[0];
rz(0.10869797) q[0];
rz(-2.3868949) q[1];
sx q[1];
rz(-1.8221868) q[1];
sx q[1];
rz(1.6356161) q[1];
rz(-0.66773141) q[2];
sx q[2];
rz(-2.3420391) q[2];
sx q[2];
rz(2.5426368) q[2];
rz(0.084372088) q[3];
sx q[3];
rz(-1.2481239) q[3];
sx q[3];
rz(0.39831755) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];