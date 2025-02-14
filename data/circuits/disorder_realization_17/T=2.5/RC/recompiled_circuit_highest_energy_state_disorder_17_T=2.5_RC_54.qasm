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
rz(-2.1276346) q[0];
sx q[0];
rz(-1.4528217) q[0];
sx q[0];
rz(-2.5435574) q[0];
rz(1.1818089) q[1];
sx q[1];
rz(-0.66507566) q[1];
sx q[1];
rz(-2.9887078) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11422296) q[0];
sx q[0];
rz(-1.682617) q[0];
sx q[0];
rz(-3.0614475) q[0];
rz(-pi) q[1];
rz(0.48798765) q[2];
sx q[2];
rz(-1.5190912) q[2];
sx q[2];
rz(-1.082513) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8522332) q[1];
sx q[1];
rz(-2.8948862) q[1];
sx q[1];
rz(-2.1092595) q[1];
rz(-0.064611994) q[3];
sx q[3];
rz(-0.046053208) q[3];
sx q[3];
rz(2.8833431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.87024706) q[2];
sx q[2];
rz(-0.18834867) q[2];
sx q[2];
rz(-1.9625473) q[2];
rz(0.41006586) q[3];
sx q[3];
rz(-1.054801) q[3];
sx q[3];
rz(0.89455354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1882741) q[0];
sx q[0];
rz(-2.4725547) q[0];
sx q[0];
rz(-2.8642995) q[0];
rz(-2.9412728) q[1];
sx q[1];
rz(-2.2453313) q[1];
sx q[1];
rz(3.0576113) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1071716) q[0];
sx q[0];
rz(-2.0382529) q[0];
sx q[0];
rz(-2.9185027) q[0];
x q[1];
rz(-1.7189156) q[2];
sx q[2];
rz(-0.84858719) q[2];
sx q[2];
rz(2.2715457) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18296227) q[1];
sx q[1];
rz(-2.0271027) q[1];
sx q[1];
rz(-2.0668304) q[1];
rz(-pi) q[2];
rz(1.8324319) q[3];
sx q[3];
rz(-2.0896834) q[3];
sx q[3];
rz(-0.46552697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1818485) q[2];
sx q[2];
rz(-2.4740348) q[2];
sx q[2];
rz(-0.7461156) q[2];
rz(2.5203846) q[3];
sx q[3];
rz(-2.4977903) q[3];
sx q[3];
rz(1.3889036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8212432) q[0];
sx q[0];
rz(-0.72750434) q[0];
sx q[0];
rz(1.2028836) q[0];
rz(2.7299643) q[1];
sx q[1];
rz(-1.2806712) q[1];
sx q[1];
rz(2.0278377) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4899556) q[0];
sx q[0];
rz(-0.86511602) q[0];
sx q[0];
rz(-0.56244992) q[0];
rz(-pi) q[1];
x q[1];
rz(2.162777) q[2];
sx q[2];
rz(-2.4360058) q[2];
sx q[2];
rz(1.5075114) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.74677709) q[1];
sx q[1];
rz(-1.703843) q[1];
sx q[1];
rz(-2.2486671) q[1];
x q[2];
rz(-0.87388523) q[3];
sx q[3];
rz(-2.0351119) q[3];
sx q[3];
rz(3.0261542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0456298) q[2];
sx q[2];
rz(-2.3489504) q[2];
sx q[2];
rz(1.9795798) q[2];
rz(2.3368733) q[3];
sx q[3];
rz(-1.8685124) q[3];
sx q[3];
rz(1.2812251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3941037) q[0];
sx q[0];
rz(-1.2825613) q[0];
sx q[0];
rz(1.5149186) q[0];
rz(-0.50846848) q[1];
sx q[1];
rz(-0.6503121) q[1];
sx q[1];
rz(-1.633684) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2861904) q[0];
sx q[0];
rz(-1.5729181) q[0];
sx q[0];
rz(-1.7218263) q[0];
rz(0.99782439) q[2];
sx q[2];
rz(-2.8854239) q[2];
sx q[2];
rz(-2.1081971) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0421988) q[1];
sx q[1];
rz(-1.6833937) q[1];
sx q[1];
rz(-1.2641175) q[1];
rz(2.1056469) q[3];
sx q[3];
rz(-0.7300762) q[3];
sx q[3];
rz(1.7022145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.41370979) q[2];
sx q[2];
rz(-1.1943613) q[2];
sx q[2];
rz(-0.74753648) q[2];
rz(-1.2306635) q[3];
sx q[3];
rz(-2.3321798) q[3];
sx q[3];
rz(2.6443853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43664765) q[0];
sx q[0];
rz(-1.1282938) q[0];
sx q[0];
rz(-1.6823912) q[0];
rz(-0.35405007) q[1];
sx q[1];
rz(-2.470128) q[1];
sx q[1];
rz(0.57463247) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45978776) q[0];
sx q[0];
rz(-1.5842251) q[0];
sx q[0];
rz(3.0643377) q[0];
rz(-pi) q[1];
rz(-0.76894297) q[2];
sx q[2];
rz(-1.3731602) q[2];
sx q[2];
rz(2.1773424) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1916442) q[1];
sx q[1];
rz(-0.23298082) q[1];
sx q[1];
rz(-0.22535546) q[1];
x q[2];
rz(-2.0784723) q[3];
sx q[3];
rz(-0.23323828) q[3];
sx q[3];
rz(2.4539029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6358801) q[2];
sx q[2];
rz(-1.1735703) q[2];
sx q[2];
rz(-0.20225784) q[2];
rz(-2.108719) q[3];
sx q[3];
rz(-2.5346916) q[3];
sx q[3];
rz(0.066298299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.031484) q[0];
sx q[0];
rz(-0.38723543) q[0];
sx q[0];
rz(-2.7701344) q[0];
rz(-0.29608852) q[1];
sx q[1];
rz(-1.5874054) q[1];
sx q[1];
rz(2.9782226) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0609157) q[0];
sx q[0];
rz(-2.0656343) q[0];
sx q[0];
rz(1.3617152) q[0];
rz(2.4642015) q[2];
sx q[2];
rz(-2.4755726) q[2];
sx q[2];
rz(0.99413727) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0381546) q[1];
sx q[1];
rz(-2.0088327) q[1];
sx q[1];
rz(-1.7704727) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7512239) q[3];
sx q[3];
rz(-0.70367763) q[3];
sx q[3];
rz(-2.0247519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0801733) q[2];
sx q[2];
rz(-0.27013186) q[2];
sx q[2];
rz(-2.488625) q[2];
rz(2.669892) q[3];
sx q[3];
rz(-2.3931849) q[3];
sx q[3];
rz(-2.4145943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2111874) q[0];
sx q[0];
rz(-2.0624332) q[0];
sx q[0];
rz(1.0373254) q[0];
rz(-2.6264181) q[1];
sx q[1];
rz(-1.8637135) q[1];
sx q[1];
rz(-1.3264664) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4600996) q[0];
sx q[0];
rz(-1.9578448) q[0];
sx q[0];
rz(-0.99868628) q[0];
rz(-pi) q[1];
rz(3.0318119) q[2];
sx q[2];
rz(-2.7876543) q[2];
sx q[2];
rz(-0.76895163) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4660037) q[1];
sx q[1];
rz(-1.6674622) q[1];
sx q[1];
rz(-0.76670353) q[1];
rz(-pi) q[2];
rz(-1.8765728) q[3];
sx q[3];
rz(-1.5103087) q[3];
sx q[3];
rz(0.6636338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4482164) q[2];
sx q[2];
rz(-2.6770112) q[2];
sx q[2];
rz(-2.8153937) q[2];
rz(2.163573) q[3];
sx q[3];
rz(-1.2468485) q[3];
sx q[3];
rz(3.0514362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.4433032) q[0];
sx q[0];
rz(-2.8841618) q[0];
sx q[0];
rz(0.28939104) q[0];
rz(-0.012271317) q[1];
sx q[1];
rz(-2.8616276) q[1];
sx q[1];
rz(1.4853005) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0146831) q[0];
sx q[0];
rz(-2.2193847) q[0];
sx q[0];
rz(0.65833433) q[0];
rz(-pi) q[1];
rz(2.907862) q[2];
sx q[2];
rz(-2.9119891) q[2];
sx q[2];
rz(0.52403852) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.030823) q[1];
sx q[1];
rz(-1.4258037) q[1];
sx q[1];
rz(0.20408634) q[1];
x q[2];
rz(-0.31987736) q[3];
sx q[3];
rz(-1.1231224) q[3];
sx q[3];
rz(0.059739865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5926008) q[2];
sx q[2];
rz(-1.5718549) q[2];
sx q[2];
rz(2.9686417) q[2];
rz(1.7671827) q[3];
sx q[3];
rz(-1.7632615) q[3];
sx q[3];
rz(-1.6636728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3953111) q[0];
sx q[0];
rz(-0.44773856) q[0];
sx q[0];
rz(-0.81004274) q[0];
rz(-2.7250302) q[1];
sx q[1];
rz(-2.3655393) q[1];
sx q[1];
rz(-0.3449482) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.024558) q[0];
sx q[0];
rz(-0.86625615) q[0];
sx q[0];
rz(-1.7323094) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0428997) q[2];
sx q[2];
rz(-2.8213892) q[2];
sx q[2];
rz(0.36123438) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7899985) q[1];
sx q[1];
rz(-1.8934544) q[1];
sx q[1];
rz(-2.0568454) q[1];
x q[2];
rz(0.79699253) q[3];
sx q[3];
rz(-0.16599645) q[3];
sx q[3];
rz(-0.57193631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.714146) q[2];
sx q[2];
rz(-0.88627187) q[2];
sx q[2];
rz(-3.0575338) q[2];
rz(-2.2345624) q[3];
sx q[3];
rz(-1.7232938) q[3];
sx q[3];
rz(2.1840054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8091938) q[0];
sx q[0];
rz(-0.087662307) q[0];
sx q[0];
rz(-2.8293389) q[0];
rz(-2.9088083) q[1];
sx q[1];
rz(-0.39118958) q[1];
sx q[1];
rz(-1.2871453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2904743) q[0];
sx q[0];
rz(-2.2016207) q[0];
sx q[0];
rz(1.9734971) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.078545) q[2];
sx q[2];
rz(-0.64008605) q[2];
sx q[2];
rz(-0.62713059) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6617468) q[1];
sx q[1];
rz(-1.8150618) q[1];
sx q[1];
rz(2.0354009) q[1];
rz(-2.5165043) q[3];
sx q[3];
rz(-1.2660789) q[3];
sx q[3];
rz(-2.8947322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1708258) q[2];
sx q[2];
rz(-1.2831186) q[2];
sx q[2];
rz(-1.6440561) q[2];
rz(1.0981285) q[3];
sx q[3];
rz(-0.69881717) q[3];
sx q[3];
rz(0.48925492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.004414) q[0];
sx q[0];
rz(-1.7902086) q[0];
sx q[0];
rz(-0.93850346) q[0];
rz(1.3068403) q[1];
sx q[1];
rz(-2.2521781) q[1];
sx q[1];
rz(2.3440012) q[1];
rz(-0.88181344) q[2];
sx q[2];
rz(-1.4588933) q[2];
sx q[2];
rz(2.3068538) q[2];
rz(0.14214235) q[3];
sx q[3];
rz(-0.90928838) q[3];
sx q[3];
rz(1.7293775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
