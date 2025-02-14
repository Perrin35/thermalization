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
rz(-1.9189605) q[0];
sx q[0];
rz(-0.55813342) q[0];
sx q[0];
rz(0.65879917) q[0];
rz(-2.2539723) q[1];
sx q[1];
rz(-2.4467111) q[1];
sx q[1];
rz(-0.08858362) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50669955) q[0];
sx q[0];
rz(-1.1453712) q[0];
sx q[0];
rz(-2.9296207) q[0];
x q[1];
rz(-2.7062662) q[2];
sx q[2];
rz(-2.0930039) q[2];
sx q[2];
rz(-0.84964035) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3693847) q[1];
sx q[1];
rz(-1.4190892) q[1];
sx q[1];
rz(0.94767344) q[1];
rz(-pi) q[2];
rz(2.3490099) q[3];
sx q[3];
rz(-1.152609) q[3];
sx q[3];
rz(-1.7824695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6170071) q[2];
sx q[2];
rz(-1.1716537) q[2];
sx q[2];
rz(0.48801547) q[2];
rz(-2.6739142) q[3];
sx q[3];
rz(-0.52507639) q[3];
sx q[3];
rz(0.16744965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2117598) q[0];
sx q[0];
rz(-0.81120315) q[0];
sx q[0];
rz(1.8875341) q[0];
rz(0.60011855) q[1];
sx q[1];
rz(-1.3328726) q[1];
sx q[1];
rz(2.7235203) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1195898) q[0];
sx q[0];
rz(-1.3922577) q[0];
sx q[0];
rz(-2.5994632) q[0];
x q[1];
rz(-0.64855021) q[2];
sx q[2];
rz(-1.6706711) q[2];
sx q[2];
rz(1.5431736) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4198735) q[1];
sx q[1];
rz(-1.5054107) q[1];
sx q[1];
rz(1.2544778) q[1];
rz(-pi) q[2];
rz(-2.2369305) q[3];
sx q[3];
rz(-2.2719529) q[3];
sx q[3];
rz(1.9405492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7347001) q[2];
sx q[2];
rz(-2.0105346) q[2];
sx q[2];
rz(-0.1184173) q[2];
rz(-2.7393869) q[3];
sx q[3];
rz(-1.6536568) q[3];
sx q[3];
rz(-1.8124628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44620946) q[0];
sx q[0];
rz(-2.8948247) q[0];
sx q[0];
rz(1.9870019) q[0];
rz(2.0788976) q[1];
sx q[1];
rz(-2.1176391) q[1];
sx q[1];
rz(-1.9564995) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0010312) q[0];
sx q[0];
rz(-0.90406448) q[0];
sx q[0];
rz(2.786955) q[0];
rz(-pi) q[1];
rz(0.22151557) q[2];
sx q[2];
rz(-2.0755526) q[2];
sx q[2];
rz(-1.5496448) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9251057) q[1];
sx q[1];
rz(-0.47678927) q[1];
sx q[1];
rz(-2.2170223) q[1];
x q[2];
rz(-0.40374741) q[3];
sx q[3];
rz(-1.8093411) q[3];
sx q[3];
rz(2.8806137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2350754) q[2];
sx q[2];
rz(-2.2791028) q[2];
sx q[2];
rz(-1.4441351) q[2];
rz(-2.2200572) q[3];
sx q[3];
rz(-0.60465616) q[3];
sx q[3];
rz(0.055189341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4398572) q[0];
sx q[0];
rz(-3.0210962) q[0];
sx q[0];
rz(-3.0635656) q[0];
rz(-2.0241418) q[1];
sx q[1];
rz(-1.8945339) q[1];
sx q[1];
rz(1.4720346) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0271281) q[0];
sx q[0];
rz(-2.3660866) q[0];
sx q[0];
rz(-1.7239611) q[0];
rz(0.83547893) q[2];
sx q[2];
rz(-2.4148921) q[2];
sx q[2];
rz(-1.5269321) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7546447) q[1];
sx q[1];
rz(-2.7803763) q[1];
sx q[1];
rz(1.4903461) q[1];
rz(-1.941664) q[3];
sx q[3];
rz(-2.2590593) q[3];
sx q[3];
rz(2.1394068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3052519) q[2];
sx q[2];
rz(-0.64771104) q[2];
sx q[2];
rz(0.15023896) q[2];
rz(-2.5495106) q[3];
sx q[3];
rz(-1.3034857) q[3];
sx q[3];
rz(1.1881812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-2.2109461) q[0];
sx q[0];
rz(-0.35363126) q[0];
sx q[0];
rz(-0.76552248) q[0];
rz(-0.80134478) q[1];
sx q[1];
rz(-1.067602) q[1];
sx q[1];
rz(1.9498922) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36574349) q[0];
sx q[0];
rz(-1.8738215) q[0];
sx q[0];
rz(-0.55935212) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4214864) q[2];
sx q[2];
rz(-0.65944797) q[2];
sx q[2];
rz(2.5044005) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0914407) q[1];
sx q[1];
rz(-1.5066083) q[1];
sx q[1];
rz(-2.2540171) q[1];
rz(-pi) q[2];
rz(1.3016316) q[3];
sx q[3];
rz(-1.3084305) q[3];
sx q[3];
rz(3.0848665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2661065) q[2];
sx q[2];
rz(-2.2849639) q[2];
sx q[2];
rz(-2.6066656) q[2];
rz(-2.5180425) q[3];
sx q[3];
rz(-2.0303191) q[3];
sx q[3];
rz(-0.88304869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9542338) q[0];
sx q[0];
rz(-2.7281902) q[0];
sx q[0];
rz(-2.2075388) q[0];
rz(1.8961228) q[1];
sx q[1];
rz(-1.5555236) q[1];
sx q[1];
rz(0.62560558) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0348822) q[0];
sx q[0];
rz(-0.83256522) q[0];
sx q[0];
rz(0.24068479) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.093354524) q[2];
sx q[2];
rz(-0.90269222) q[2];
sx q[2];
rz(-2.7808288) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.44613353) q[1];
sx q[1];
rz(-1.2810105) q[1];
sx q[1];
rz(-0.70420806) q[1];
rz(-pi) q[2];
rz(-0.91208338) q[3];
sx q[3];
rz(-0.28277031) q[3];
sx q[3];
rz(1.2706437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.049456747) q[2];
sx q[2];
rz(-0.9149887) q[2];
sx q[2];
rz(0.12506872) q[2];
rz(-0.72758979) q[3];
sx q[3];
rz(-1.5743419) q[3];
sx q[3];
rz(-2.6653813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8759988) q[0];
sx q[0];
rz(-2.6662874) q[0];
sx q[0];
rz(1.884961) q[0];
rz(1.1988634) q[1];
sx q[1];
rz(-1.4224982) q[1];
sx q[1];
rz(1.2748324) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.066028) q[0];
sx q[0];
rz(-2.4504553) q[0];
sx q[0];
rz(-0.71623556) q[0];
rz(1.0339678) q[2];
sx q[2];
rz(-1.8933927) q[2];
sx q[2];
rz(-0.19746298) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.82060888) q[1];
sx q[1];
rz(-0.57772355) q[1];
sx q[1];
rz(2.251365) q[1];
rz(1.2009462) q[3];
sx q[3];
rz(-0.33888926) q[3];
sx q[3];
rz(1.8367934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4510497) q[2];
sx q[2];
rz(-1.0131016) q[2];
sx q[2];
rz(0.96735442) q[2];
rz(1.5830154) q[3];
sx q[3];
rz(-1.5019006) q[3];
sx q[3];
rz(-0.30174747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79307443) q[0];
sx q[0];
rz(-0.59978849) q[0];
sx q[0];
rz(-3.0317958) q[0];
rz(-1.173136) q[1];
sx q[1];
rz(-0.99183142) q[1];
sx q[1];
rz(-2.0593624) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30352509) q[0];
sx q[0];
rz(-2.3278624) q[0];
sx q[0];
rz(2.8007048) q[0];
x q[1];
rz(2.0556446) q[2];
sx q[2];
rz(-1.4746801) q[2];
sx q[2];
rz(2.3342867) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6761192) q[1];
sx q[1];
rz(-2.2567368) q[1];
sx q[1];
rz(0.53829389) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9475329) q[3];
sx q[3];
rz(-2.1749745) q[3];
sx q[3];
rz(-2.8967711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.75862306) q[2];
sx q[2];
rz(-2.0489645) q[2];
sx q[2];
rz(0.36337241) q[2];
rz(0.97909561) q[3];
sx q[3];
rz(-2.4978814) q[3];
sx q[3];
rz(2.6399829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30524224) q[0];
sx q[0];
rz(-0.79638052) q[0];
sx q[0];
rz(-0.35261944) q[0];
rz(1.3615707) q[1];
sx q[1];
rz(-1.1905328) q[1];
sx q[1];
rz(-3.000066) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82194009) q[0];
sx q[0];
rz(-1.7403649) q[0];
sx q[0];
rz(1.5362306) q[0];
x q[1];
rz(2.5472801) q[2];
sx q[2];
rz(-1.9611352) q[2];
sx q[2];
rz(0.28536404) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80672164) q[1];
sx q[1];
rz(-1.9892684) q[1];
sx q[1];
rz(-1.970402) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0120506) q[3];
sx q[3];
rz(-2.1390474) q[3];
sx q[3];
rz(0.38145867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8286459) q[2];
sx q[2];
rz(-1.6673648) q[2];
sx q[2];
rz(-2.055577) q[2];
rz(-2.9588251) q[3];
sx q[3];
rz(-2.1153617) q[3];
sx q[3];
rz(2.1695547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52520853) q[0];
sx q[0];
rz(-1.5807736) q[0];
sx q[0];
rz(-0.64424789) q[0];
rz(-3.0746025) q[1];
sx q[1];
rz(-1.3863775) q[1];
sx q[1];
rz(3.0279874) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85746914) q[0];
sx q[0];
rz(-2.9175287) q[0];
sx q[0];
rz(1.5518918) q[0];
rz(-pi) q[1];
rz(-0.082661672) q[2];
sx q[2];
rz(-1.6726521) q[2];
sx q[2];
rz(1.5422623) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.41109243) q[1];
sx q[1];
rz(-2.9343795) q[1];
sx q[1];
rz(-2.5403008) q[1];
x q[2];
rz(-1.4573649) q[3];
sx q[3];
rz(-0.3598752) q[3];
sx q[3];
rz(-2.7168093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.99864787) q[2];
sx q[2];
rz(-1.8871658) q[2];
sx q[2];
rz(0.76510731) q[2];
rz(0.1782002) q[3];
sx q[3];
rz(-1.5604115) q[3];
sx q[3];
rz(0.83711886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9255623) q[0];
sx q[0];
rz(-2.3860274) q[0];
sx q[0];
rz(-0.26269333) q[0];
rz(-0.33949159) q[1];
sx q[1];
rz(-1.2295634) q[1];
sx q[1];
rz(-2.0801574) q[1];
rz(-0.65503623) q[2];
sx q[2];
rz(-1.8476386) q[2];
sx q[2];
rz(2.4551433) q[2];
rz(2.8398561) q[3];
sx q[3];
rz(-1.269125) q[3];
sx q[3];
rz(2.453809) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
