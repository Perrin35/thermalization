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
rz(-1.2008774) q[0];
sx q[0];
rz(-0.90056363) q[0];
sx q[0];
rz(-0.23298921) q[0];
rz(-1.5268582) q[1];
sx q[1];
rz(-2.0695217) q[1];
sx q[1];
rz(1.1195247) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096123556) q[0];
sx q[0];
rz(-1.5579559) q[0];
sx q[0];
rz(-1.8052285) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.018466516) q[2];
sx q[2];
rz(-2.5800843) q[2];
sx q[2];
rz(2.7811921) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.865316) q[1];
sx q[1];
rz(-2.7807693) q[1];
sx q[1];
rz(-1.7860543) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5324838) q[3];
sx q[3];
rz(-1.850633) q[3];
sx q[3];
rz(2.5697054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4423674) q[2];
sx q[2];
rz(-1.0590326) q[2];
sx q[2];
rz(-3.0287058) q[2];
rz(2.8804307) q[3];
sx q[3];
rz(-1.7966725) q[3];
sx q[3];
rz(-1.8908709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79171044) q[0];
sx q[0];
rz(-3.0630906) q[0];
sx q[0];
rz(0.054340266) q[0];
rz(0.21866523) q[1];
sx q[1];
rz(-1.668914) q[1];
sx q[1];
rz(2.7770619) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69037777) q[0];
sx q[0];
rz(-1.0752819) q[0];
sx q[0];
rz(2.4507603) q[0];
rz(-pi) q[1];
rz(-0.40016115) q[2];
sx q[2];
rz(-0.98099835) q[2];
sx q[2];
rz(-1.9831744) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6333225) q[1];
sx q[1];
rz(-0.62707096) q[1];
sx q[1];
rz(3.0116664) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9596635) q[3];
sx q[3];
rz(-1.9046138) q[3];
sx q[3];
rz(0.16063375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.24078044) q[2];
sx q[2];
rz(-1.4419) q[2];
sx q[2];
rz(-2.0330632) q[2];
rz(2.945914) q[3];
sx q[3];
rz(-1.207573) q[3];
sx q[3];
rz(1.4669363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7473258) q[0];
sx q[0];
rz(-1.0402004) q[0];
sx q[0];
rz(-2.105383) q[0];
rz(-0.58468435) q[1];
sx q[1];
rz(-1.4671289) q[1];
sx q[1];
rz(-0.55589693) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5591878) q[0];
sx q[0];
rz(-0.83600658) q[0];
sx q[0];
rz(0.57484267) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1039882) q[2];
sx q[2];
rz(-1.4935263) q[2];
sx q[2];
rz(-1.8118878) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8396254) q[1];
sx q[1];
rz(-2.4263067) q[1];
sx q[1];
rz(1.3705181) q[1];
x q[2];
rz(0.44938748) q[3];
sx q[3];
rz(-1.3419749) q[3];
sx q[3];
rz(0.23923161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3802152) q[2];
sx q[2];
rz(-0.20647241) q[2];
sx q[2];
rz(1.9492487) q[2];
rz(-0.10382593) q[3];
sx q[3];
rz(-1.1622279) q[3];
sx q[3];
rz(1.9942358) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1059859) q[0];
sx q[0];
rz(-0.62653956) q[0];
sx q[0];
rz(-1.71738) q[0];
rz(-2.4183938) q[1];
sx q[1];
rz(-2.9056748) q[1];
sx q[1];
rz(-2.9211488) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0766692) q[0];
sx q[0];
rz(-0.60416302) q[0];
sx q[0];
rz(2.2518065) q[0];
rz(2.0218684) q[2];
sx q[2];
rz(-2.5325518) q[2];
sx q[2];
rz(0.90897564) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.087763) q[1];
sx q[1];
rz(-1.9265895) q[1];
sx q[1];
rz(0.1202113) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1573767) q[3];
sx q[3];
rz(-1.2261021) q[3];
sx q[3];
rz(-0.34806309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.38349884) q[2];
sx q[2];
rz(-1.8269962) q[2];
sx q[2];
rz(0.51551762) q[2];
rz(-0.085144194) q[3];
sx q[3];
rz(-0.65535039) q[3];
sx q[3];
rz(2.3084124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4715217) q[0];
sx q[0];
rz(-0.18297289) q[0];
sx q[0];
rz(-2.4772189) q[0];
rz(-2.6145256) q[1];
sx q[1];
rz(-1.138569) q[1];
sx q[1];
rz(1.9648431) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.949958) q[0];
sx q[0];
rz(-1.598888) q[0];
sx q[0];
rz(-1.4856443) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46026261) q[2];
sx q[2];
rz(-2.2027594) q[2];
sx q[2];
rz(3.021701) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0461573) q[1];
sx q[1];
rz(-1.5514138) q[1];
sx q[1];
rz(-2.9033015) q[1];
rz(-pi) q[2];
rz(1.845729) q[3];
sx q[3];
rz(-0.58919135) q[3];
sx q[3];
rz(-1.2146035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.027017) q[2];
sx q[2];
rz(-0.21275529) q[2];
sx q[2];
rz(0.93811718) q[2];
rz(-2.0959057) q[3];
sx q[3];
rz(-0.18200471) q[3];
sx q[3];
rz(-0.81030455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.806458) q[0];
sx q[0];
rz(-1.198575) q[0];
sx q[0];
rz(1.2499811) q[0];
rz(-0.20420034) q[1];
sx q[1];
rz(-2.4649492) q[1];
sx q[1];
rz(-2.2883889) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9850904) q[0];
sx q[0];
rz(-1.2285875) q[0];
sx q[0];
rz(-1.2669105) q[0];
rz(-0.95909024) q[2];
sx q[2];
rz(-1.0420799) q[2];
sx q[2];
rz(2.6564646) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2320441) q[1];
sx q[1];
rz(-2.2799788) q[1];
sx q[1];
rz(-2.1348165) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9123934) q[3];
sx q[3];
rz(-1.4506243) q[3];
sx q[3];
rz(-0.12334331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.044067232) q[2];
sx q[2];
rz(-1.7996457) q[2];
sx q[2];
rz(1.9419144) q[2];
rz(-2.1357644) q[3];
sx q[3];
rz(-2.2483716) q[3];
sx q[3];
rz(-0.83542663) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47578874) q[0];
sx q[0];
rz(-0.26837334) q[0];
sx q[0];
rz(-2.1589808) q[0];
rz(1.3081374) q[1];
sx q[1];
rz(-1.9255226) q[1];
sx q[1];
rz(-1.7024202) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84837259) q[0];
sx q[0];
rz(-1.2407899) q[0];
sx q[0];
rz(0.35949913) q[0];
rz(-pi) q[1];
rz(2.0928395) q[2];
sx q[2];
rz(-1.6380651) q[2];
sx q[2];
rz(-1.0053588) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8285259) q[1];
sx q[1];
rz(-1.8983574) q[1];
sx q[1];
rz(0.86159535) q[1];
rz(-1.7129219) q[3];
sx q[3];
rz(-2.0439337) q[3];
sx q[3];
rz(2.9700235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9787489) q[2];
sx q[2];
rz(-2.2480201) q[2];
sx q[2];
rz(2.5992375) q[2];
rz(2.1584623) q[3];
sx q[3];
rz(-1.977908) q[3];
sx q[3];
rz(-2.2014309) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56383175) q[0];
sx q[0];
rz(-1.4262119) q[0];
sx q[0];
rz(-2.3610709) q[0];
rz(1.966194) q[1];
sx q[1];
rz(-0.54882097) q[1];
sx q[1];
rz(1.0343879) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7808409) q[0];
sx q[0];
rz(-0.409161) q[0];
sx q[0];
rz(0.43171127) q[0];
x q[1];
rz(-1.0532736) q[2];
sx q[2];
rz(-2.9773601) q[2];
sx q[2];
rz(2.9688778) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3745034) q[1];
sx q[1];
rz(-1.6853635) q[1];
sx q[1];
rz(0.40278816) q[1];
rz(-pi) q[2];
rz(-0.42203776) q[3];
sx q[3];
rz(-2.5521818) q[3];
sx q[3];
rz(0.60822076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9913651) q[2];
sx q[2];
rz(-1.4069724) q[2];
sx q[2];
rz(-0.11037174) q[2];
rz(1.8433833) q[3];
sx q[3];
rz(-1.7896264) q[3];
sx q[3];
rz(-2.8492294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2503535) q[0];
sx q[0];
rz(-2.1284916) q[0];
sx q[0];
rz(1.891834) q[0];
rz(-3.0505782) q[1];
sx q[1];
rz(-2.4200771) q[1];
sx q[1];
rz(-2.4299842) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7994493) q[0];
sx q[0];
rz(-0.79252386) q[0];
sx q[0];
rz(-1.5802797) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3468412) q[2];
sx q[2];
rz(-2.134179) q[2];
sx q[2];
rz(1.7573485) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6775637) q[1];
sx q[1];
rz(-1.11598) q[1];
sx q[1];
rz(1.1879261) q[1];
rz(-1.2290088) q[3];
sx q[3];
rz(-0.84484378) q[3];
sx q[3];
rz(-1.0085229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.3332425) q[2];
sx q[2];
rz(-1.658172) q[2];
sx q[2];
rz(-1.9688155) q[2];
rz(-1.5277537) q[3];
sx q[3];
rz(-1.0781735) q[3];
sx q[3];
rz(-0.095002256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.461819) q[0];
sx q[0];
rz(-3.0175896) q[0];
sx q[0];
rz(-1.5754196) q[0];
rz(-0.35789403) q[1];
sx q[1];
rz(-2.0707668) q[1];
sx q[1];
rz(-0.90248743) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6704412) q[0];
sx q[0];
rz(-0.67865279) q[0];
sx q[0];
rz(-1.9891692) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4975791) q[2];
sx q[2];
rz(-1.5884382) q[2];
sx q[2];
rz(-2.0669075) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5229228) q[1];
sx q[1];
rz(-2.3176207) q[1];
sx q[1];
rz(2.0342779) q[1];
rz(1.020458) q[3];
sx q[3];
rz(-1.6452879) q[3];
sx q[3];
rz(1.0005152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1050528) q[2];
sx q[2];
rz(-2.5532711) q[2];
sx q[2];
rz(-1.499929) q[2];
rz(-0.52122742) q[3];
sx q[3];
rz(-0.14385496) q[3];
sx q[3];
rz(2.1874793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70984107) q[0];
sx q[0];
rz(-2.3713645) q[0];
sx q[0];
rz(1.4159528) q[0];
rz(-0.00090986666) q[1];
sx q[1];
rz(-1.6699427) q[1];
sx q[1];
rz(-1.4572399) q[1];
rz(-1.7224689) q[2];
sx q[2];
rz(-2.5628452) q[2];
sx q[2];
rz(2.645523) q[2];
rz(0.88738995) q[3];
sx q[3];
rz(-0.94684623) q[3];
sx q[3];
rz(2.0317247) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
