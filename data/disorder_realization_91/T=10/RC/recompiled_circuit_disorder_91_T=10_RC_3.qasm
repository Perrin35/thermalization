OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.70541731) q[0];
sx q[0];
rz(-2.5751312) q[0];
sx q[0];
rz(-0.17106549) q[0];
rz(1.45362) q[1];
sx q[1];
rz(-0.34314081) q[1];
sx q[1];
rz(-1.3309825) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1216461) q[0];
sx q[0];
rz(-2.2197147) q[0];
sx q[0];
rz(-1.5050423) q[0];
rz(-pi) q[1];
rz(1.5924442) q[2];
sx q[2];
rz(-1.0550371) q[2];
sx q[2];
rz(1.5390918) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0806182) q[1];
sx q[1];
rz(-1.8218826) q[1];
sx q[1];
rz(-3.1239448) q[1];
x q[2];
rz(2.8475548) q[3];
sx q[3];
rz(-0.66411823) q[3];
sx q[3];
rz(0.77120632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3657637) q[2];
sx q[2];
rz(-0.77314955) q[2];
sx q[2];
rz(1.8120871) q[2];
rz(1.7261516) q[3];
sx q[3];
rz(-1.5664145) q[3];
sx q[3];
rz(-2.9392488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1502894) q[0];
sx q[0];
rz(-0.54420272) q[0];
sx q[0];
rz(-1.4527028) q[0];
rz(2.9406722) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(2.8413049) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072104134) q[0];
sx q[0];
rz(-2.527643) q[0];
sx q[0];
rz(2.9314562) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57057256) q[2];
sx q[2];
rz(-1.3010446) q[2];
sx q[2];
rz(-0.74769339) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2806432) q[1];
sx q[1];
rz(-1.7935392) q[1];
sx q[1];
rz(-1.8722948) q[1];
rz(1.7598011) q[3];
sx q[3];
rz(-1.731589) q[3];
sx q[3];
rz(0.27634987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.77256569) q[2];
sx q[2];
rz(-0.3521266) q[2];
sx q[2];
rz(-1.0380113) q[2];
rz(-2.299262) q[3];
sx q[3];
rz(-1.7269644) q[3];
sx q[3];
rz(0.37030181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5623986) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(2.753479) q[0];
rz(0.072470486) q[1];
sx q[1];
rz(-1.4265172) q[1];
sx q[1];
rz(0.31633502) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4916723) q[0];
sx q[0];
rz(-0.24501093) q[0];
sx q[0];
rz(-1.578376) q[0];
rz(-pi) q[1];
rz(-0.6109654) q[2];
sx q[2];
rz(-0.66662153) q[2];
sx q[2];
rz(-1.1952458) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0132732) q[1];
sx q[1];
rz(-1.1228021) q[1];
sx q[1];
rz(-1.1303348) q[1];
rz(-pi) q[2];
rz(-1.2315138) q[3];
sx q[3];
rz(-1.2959891) q[3];
sx q[3];
rz(-2.2846082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0324273) q[2];
sx q[2];
rz(-2.9563603) q[2];
sx q[2];
rz(-0.16080984) q[2];
rz(-0.12604776) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(-1.0236615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2407103) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(-2.3098992) q[0];
rz(1.7968934) q[1];
sx q[1];
rz(-2.1422377) q[1];
sx q[1];
rz(1.6279189) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9261525) q[0];
sx q[0];
rz(-2.2427796) q[0];
sx q[0];
rz(-1.6295022) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0604826) q[2];
sx q[2];
rz(-2.1334279) q[2];
sx q[2];
rz(0.62603355) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1308243) q[1];
sx q[1];
rz(-2.0190034) q[1];
sx q[1];
rz(-0.018647714) q[1];
rz(-2.4098445) q[3];
sx q[3];
rz(-1.2536612) q[3];
sx q[3];
rz(1.5750615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7164798) q[2];
sx q[2];
rz(-1.0093062) q[2];
sx q[2];
rz(1.9173737) q[2];
rz(2.7246357) q[3];
sx q[3];
rz(-1.0427534) q[3];
sx q[3];
rz(-0.52880374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083374627) q[0];
sx q[0];
rz(-2.871802) q[0];
sx q[0];
rz(-1.8378687) q[0];
rz(-0.45571348) q[1];
sx q[1];
rz(-0.2625176) q[1];
sx q[1];
rz(-0.051503332) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6131825) q[0];
sx q[0];
rz(-0.15325704) q[0];
sx q[0];
rz(2.2806703) q[0];
x q[1];
rz(-2.0578458) q[2];
sx q[2];
rz(-2.8205928) q[2];
sx q[2];
rz(3.0082862) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3225973) q[1];
sx q[1];
rz(-1.7261191) q[1];
sx q[1];
rz(2.0218693) q[1];
x q[2];
rz(1.6938985) q[3];
sx q[3];
rz(-0.8408635) q[3];
sx q[3];
rz(0.52809944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6872528) q[2];
sx q[2];
rz(-1.3663102) q[2];
sx q[2];
rz(2.5435737) q[2];
rz(-0.55001843) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(-2.0050744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(1.0705868) q[0];
sx q[0];
rz(-0.07659176) q[0];
sx q[0];
rz(-0.20198527) q[0];
rz(-2.1760991) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(-3.0583256) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.252713) q[0];
sx q[0];
rz(-1.3759334) q[0];
sx q[0];
rz(1.7478419) q[0];
rz(-pi) q[1];
rz(1.303057) q[2];
sx q[2];
rz(-1.3053075) q[2];
sx q[2];
rz(-0.43648411) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6492918) q[1];
sx q[1];
rz(-1.2052844) q[1];
sx q[1];
rz(0.11458061) q[1];
rz(-pi) q[2];
rz(-2.5842651) q[3];
sx q[3];
rz(-1.8263655) q[3];
sx q[3];
rz(0.42423466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.10963708) q[2];
sx q[2];
rz(-1.6836616) q[2];
sx q[2];
rz(-1.7588245) q[2];
rz(0.2494732) q[3];
sx q[3];
rz(-1.8981257) q[3];
sx q[3];
rz(1.680254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7145342) q[0];
sx q[0];
rz(-0.80474168) q[0];
sx q[0];
rz(0.89685857) q[0];
rz(-2.8915021) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(2.0239963) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5516978) q[0];
sx q[0];
rz(-0.74001827) q[0];
sx q[0];
rz(-2.9445573) q[0];
rz(-2.2824077) q[2];
sx q[2];
rz(-1.7297598) q[2];
sx q[2];
rz(2.8189916) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7270131) q[1];
sx q[1];
rz(-1.0158744) q[1];
sx q[1];
rz(-1.8068061) q[1];
rz(-pi) q[2];
rz(0.53447978) q[3];
sx q[3];
rz(-1.8090873) q[3];
sx q[3];
rz(-2.0508545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0573132) q[2];
sx q[2];
rz(-1.7417615) q[2];
sx q[2];
rz(-2.9220667) q[2];
rz(2.7548742) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(-0.99532551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052208386) q[0];
sx q[0];
rz(-1.5970255) q[0];
sx q[0];
rz(-0.21729939) q[0];
rz(-0.030933881) q[1];
sx q[1];
rz(-0.63806454) q[1];
sx q[1];
rz(2.4826179) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7227104) q[0];
sx q[0];
rz(-1.5960777) q[0];
sx q[0];
rz(0.77194571) q[0];
rz(-pi) q[1];
rz(-0.73924139) q[2];
sx q[2];
rz(-0.87202245) q[2];
sx q[2];
rz(3.0657363) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.79170376) q[1];
sx q[1];
rz(-0.83194299) q[1];
sx q[1];
rz(-1.8658584) q[1];
x q[2];
rz(-2.213845) q[3];
sx q[3];
rz(-1.5156931) q[3];
sx q[3];
rz(2.7411214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4387681) q[2];
sx q[2];
rz(-1.7610901) q[2];
sx q[2];
rz(0.32021114) q[2];
rz(1.6163588) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(-1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76072389) q[0];
sx q[0];
rz(-2.9601233) q[0];
sx q[0];
rz(0.6828126) q[0];
rz(-1.0702417) q[1];
sx q[1];
rz(-1.8990592) q[1];
sx q[1];
rz(-1.210093) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72737981) q[0];
sx q[0];
rz(-0.41398898) q[0];
sx q[0];
rz(-1.0340286) q[0];
rz(-pi) q[1];
rz(-1.9270883) q[2];
sx q[2];
rz(-0.63535832) q[2];
sx q[2];
rz(2.2941342) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.71483892) q[1];
sx q[1];
rz(-0.73298448) q[1];
sx q[1];
rz(-2.0183802) q[1];
x q[2];
rz(0.50001486) q[3];
sx q[3];
rz(-2.4588636) q[3];
sx q[3];
rz(-2.6422215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.575763) q[2];
sx q[2];
rz(-2.1255707) q[2];
sx q[2];
rz(-2.5749717) q[2];
rz(-2.2120655) q[3];
sx q[3];
rz(-2.1598787) q[3];
sx q[3];
rz(1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41314769) q[0];
sx q[0];
rz(-0.25756535) q[0];
sx q[0];
rz(-1.7096827) q[0];
rz(0.52629772) q[1];
sx q[1];
rz(-2.6142575) q[1];
sx q[1];
rz(0.79968232) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.66946) q[0];
sx q[0];
rz(-1.9181983) q[0];
sx q[0];
rz(-0.090263788) q[0];
rz(1.4893555) q[2];
sx q[2];
rz(-1.7310206) q[2];
sx q[2];
rz(-0.83934957) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77196808) q[1];
sx q[1];
rz(-2.5110285) q[1];
sx q[1];
rz(2.1715013) q[1];
rz(-0.47761376) q[3];
sx q[3];
rz(-0.16371809) q[3];
sx q[3];
rz(-1.175566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8563103) q[2];
sx q[2];
rz(-0.44390634) q[2];
sx q[2];
rz(-2.4463859) q[2];
rz(-0.55784145) q[3];
sx q[3];
rz(-1.3643684) q[3];
sx q[3];
rz(0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0556864) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(-1.862539) q[1];
sx q[1];
rz(-0.74395724) q[1];
sx q[1];
rz(-0.67768135) q[1];
rz(2.8056801) q[2];
sx q[2];
rz(-1.5716142) q[2];
sx q[2];
rz(-1.4128528) q[2];
rz(2.9813319) q[3];
sx q[3];
rz(-2.3370623) q[3];
sx q[3];
rz(-0.2094895) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
