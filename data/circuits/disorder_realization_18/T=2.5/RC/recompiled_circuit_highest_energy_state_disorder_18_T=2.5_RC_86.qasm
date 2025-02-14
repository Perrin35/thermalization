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
rz(1.2226322) q[0];
sx q[0];
rz(-2.5834592) q[0];
sx q[0];
rz(2.4827935) q[0];
rz(-2.2539723) q[1];
sx q[1];
rz(-2.4467111) q[1];
sx q[1];
rz(-0.08858362) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9889116) q[0];
sx q[0];
rz(-1.7636239) q[0];
sx q[0];
rz(1.1367984) q[0];
x q[1];
rz(-0.43532643) q[2];
sx q[2];
rz(-1.0485888) q[2];
sx q[2];
rz(-0.84964035) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0483612) q[1];
sx q[1];
rz(-2.1856896) q[1];
sx q[1];
rz(2.9555066) q[1];
rz(-pi) q[2];
rz(0.5578868) q[3];
sx q[3];
rz(-0.87432623) q[3];
sx q[3];
rz(0.59244746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52458557) q[2];
sx q[2];
rz(-1.1716537) q[2];
sx q[2];
rz(2.6535772) q[2];
rz(-0.46767849) q[3];
sx q[3];
rz(-0.52507639) q[3];
sx q[3];
rz(2.974143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2117598) q[0];
sx q[0];
rz(-0.81120315) q[0];
sx q[0];
rz(1.8875341) q[0];
rz(-2.5414741) q[1];
sx q[1];
rz(-1.80872) q[1];
sx q[1];
rz(-2.7235203) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.583823) q[0];
sx q[0];
rz(-1.038214) q[0];
sx q[0];
rz(1.7784276) q[0];
rz(0.64855021) q[2];
sx q[2];
rz(-1.4709215) q[2];
sx q[2];
rz(1.5431736) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.1723056) q[1];
sx q[1];
rz(-1.2551771) q[1];
sx q[1];
rz(0.068788485) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2369305) q[3];
sx q[3];
rz(-0.86963973) q[3];
sx q[3];
rz(1.9405492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40689251) q[2];
sx q[2];
rz(-1.1310581) q[2];
sx q[2];
rz(3.0231754) q[2];
rz(0.40220574) q[3];
sx q[3];
rz(-1.4879358) q[3];
sx q[3];
rz(-1.3291298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44620946) q[0];
sx q[0];
rz(-0.24676794) q[0];
sx q[0];
rz(1.1545908) q[0];
rz(-2.0788976) q[1];
sx q[1];
rz(-1.0239536) q[1];
sx q[1];
rz(-1.9564995) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20511928) q[0];
sx q[0];
rz(-1.2944049) q[0];
sx q[0];
rz(0.87262459) q[0];
rz(-1.1923033) q[2];
sx q[2];
rz(-0.54735294) q[2];
sx q[2];
rz(-1.1138091) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51260468) q[1];
sx q[1];
rz(-1.9459241) q[1];
sx q[1];
rz(2.8400322) q[1];
rz(-pi) q[2];
rz(2.7378452) q[3];
sx q[3];
rz(-1.8093411) q[3];
sx q[3];
rz(-0.26097894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9065173) q[2];
sx q[2];
rz(-2.2791028) q[2];
sx q[2];
rz(-1.6974576) q[2];
rz(0.92153543) q[3];
sx q[3];
rz(-0.60465616) q[3];
sx q[3];
rz(-3.0864033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4398572) q[0];
sx q[0];
rz(-0.12049645) q[0];
sx q[0];
rz(-0.078027092) q[0];
rz(-2.0241418) q[1];
sx q[1];
rz(-1.8945339) q[1];
sx q[1];
rz(-1.6695581) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11446458) q[0];
sx q[0];
rz(-0.77550602) q[0];
sx q[0];
rz(-1.7239611) q[0];
rz(2.3061137) q[2];
sx q[2];
rz(-0.72670055) q[2];
sx q[2];
rz(1.6146605) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.3009744) q[1];
sx q[1];
rz(-1.930791) q[1];
sx q[1];
rz(-3.1112413) q[1];
rz(2.4185798) q[3];
sx q[3];
rz(-1.8545056) q[3];
sx q[3];
rz(-0.32645389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3052519) q[2];
sx q[2];
rz(-2.4938816) q[2];
sx q[2];
rz(0.15023896) q[2];
rz(-0.59208208) q[3];
sx q[3];
rz(-1.838107) q[3];
sx q[3];
rz(-1.9534115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93064654) q[0];
sx q[0];
rz(-2.7879614) q[0];
sx q[0];
rz(-2.3760702) q[0];
rz(2.3402479) q[1];
sx q[1];
rz(-1.067602) q[1];
sx q[1];
rz(-1.1917005) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.649851) q[0];
sx q[0];
rz(-2.5132127) q[0];
sx q[0];
rz(0.5324441) q[0];
rz(-0.72010626) q[2];
sx q[2];
rz(-2.4821447) q[2];
sx q[2];
rz(0.63719213) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6995865) q[1];
sx q[1];
rz(-0.68574673) q[1];
sx q[1];
rz(1.4693292) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27169835) q[3];
sx q[3];
rz(-1.3110597) q[3];
sx q[3];
rz(-1.5560957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2661065) q[2];
sx q[2];
rz(-2.2849639) q[2];
sx q[2];
rz(-2.6066656) q[2];
rz(0.62355012) q[3];
sx q[3];
rz(-1.1112735) q[3];
sx q[3];
rz(-2.258544) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18735886) q[0];
sx q[0];
rz(-2.7281902) q[0];
sx q[0];
rz(-0.9340539) q[0];
rz(1.2454698) q[1];
sx q[1];
rz(-1.586069) q[1];
sx q[1];
rz(0.62560558) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0348822) q[0];
sx q[0];
rz(-2.3090274) q[0];
sx q[0];
rz(-2.9009079) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0482381) q[2];
sx q[2];
rz(-2.2389004) q[2];
sx q[2];
rz(-0.36076383) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6924507) q[1];
sx q[1];
rz(-2.3896273) q[1];
sx q[1];
rz(2.7099931) q[1];
x q[2];
rz(-0.17600893) q[3];
sx q[3];
rz(-1.3483244) q[3];
sx q[3];
rz(0.59222082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.049456747) q[2];
sx q[2];
rz(-2.226604) q[2];
sx q[2];
rz(-0.12506872) q[2];
rz(-0.72758979) q[3];
sx q[3];
rz(-1.5743419) q[3];
sx q[3];
rz(-2.6653813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8759988) q[0];
sx q[0];
rz(-2.6662874) q[0];
sx q[0];
rz(-1.884961) q[0];
rz(-1.9427293) q[1];
sx q[1];
rz(-1.7190944) q[1];
sx q[1];
rz(-1.2748324) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2293573) q[0];
sx q[0];
rz(-1.0692406) q[0];
sx q[0];
rz(-2.0683505) q[0];
rz(-pi) q[1];
rz(-2.1076249) q[2];
sx q[2];
rz(-1.8933927) q[2];
sx q[2];
rz(-0.19746298) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.15423552) q[1];
sx q[1];
rz(-1.2200095) q[1];
sx q[1];
rz(1.1018176) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.253264) q[3];
sx q[3];
rz(-1.6912563) q[3];
sx q[3];
rz(-2.5250556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4510497) q[2];
sx q[2];
rz(-1.0131016) q[2];
sx q[2];
rz(2.1742382) q[2];
rz(-1.5585772) q[3];
sx q[3];
rz(-1.6396921) q[3];
sx q[3];
rz(0.30174747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5062456) q[0];
sx q[0];
rz(-1.8162586) q[0];
sx q[0];
rz(2.3574791) q[0];
rz(-pi) q[1];
rz(-2.0556446) q[2];
sx q[2];
rz(-1.4746801) q[2];
sx q[2];
rz(-2.3342867) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6746862) q[1];
sx q[1];
rz(-1.162858) q[1];
sx q[1];
rz(2.3323184) q[1];
rz(-pi) q[2];
rz(-1.9475329) q[3];
sx q[3];
rz(-0.96661813) q[3];
sx q[3];
rz(2.8967711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.75862306) q[2];
sx q[2];
rz(-1.0926282) q[2];
sx q[2];
rz(0.36337241) q[2];
rz(-0.97909561) q[3];
sx q[3];
rz(-0.64371124) q[3];
sx q[3];
rz(2.6399829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8363504) q[0];
sx q[0];
rz(-0.79638052) q[0];
sx q[0];
rz(-2.7889732) q[0];
rz(-1.7800219) q[1];
sx q[1];
rz(-1.9510599) q[1];
sx q[1];
rz(3.000066) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74302076) q[0];
sx q[0];
rz(-1.5367265) q[0];
sx q[0];
rz(-0.16966804) q[0];
rz(-2.5472801) q[2];
sx q[2];
rz(-1.1804575) q[2];
sx q[2];
rz(-2.8562286) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80672164) q[1];
sx q[1];
rz(-1.9892684) q[1];
sx q[1];
rz(-1.970402) q[1];
rz(2.4961595) q[3];
sx q[3];
rz(-1.107599) q[3];
sx q[3];
rz(-2.2768159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3129468) q[2];
sx q[2];
rz(-1.6673648) q[2];
sx q[2];
rz(-2.055577) q[2];
rz(2.9588251) q[3];
sx q[3];
rz(-2.1153617) q[3];
sx q[3];
rz(-2.1695547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6163841) q[0];
sx q[0];
rz(-1.560819) q[0];
sx q[0];
rz(2.4973448) q[0];
rz(-3.0746025) q[1];
sx q[1];
rz(-1.7552152) q[1];
sx q[1];
rz(0.11360528) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2841235) q[0];
sx q[0];
rz(-2.9175287) q[0];
sx q[0];
rz(1.5897008) q[0];
rz(-1.6729987) q[2];
sx q[2];
rz(-1.6530286) q[2];
sx q[2];
rz(-3.1046347) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7305002) q[1];
sx q[1];
rz(-0.20721315) q[1];
sx q[1];
rz(2.5403008) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.213041) q[3];
sx q[3];
rz(-1.610667) q[3];
sx q[3];
rz(-1.8893582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.99864787) q[2];
sx q[2];
rz(-1.2544268) q[2];
sx q[2];
rz(-0.76510731) q[2];
rz(2.9633925) q[3];
sx q[3];
rz(-1.5811812) q[3];
sx q[3];
rz(-2.3044738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2160303) q[0];
sx q[0];
rz(-2.3860274) q[0];
sx q[0];
rz(-0.26269333) q[0];
rz(2.8021011) q[1];
sx q[1];
rz(-1.2295634) q[1];
sx q[1];
rz(-2.0801574) q[1];
rz(-2.4865564) q[2];
sx q[2];
rz(-1.2939541) q[2];
sx q[2];
rz(-0.68644938) q[2];
rz(0.30173652) q[3];
sx q[3];
rz(-1.8724676) q[3];
sx q[3];
rz(-0.68778366) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
