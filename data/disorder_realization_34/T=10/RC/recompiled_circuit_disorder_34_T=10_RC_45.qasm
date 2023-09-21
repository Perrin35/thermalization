OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.37378398) q[0];
sx q[0];
rz(-2.7019579) q[0];
sx q[0];
rz(0.08131942) q[0];
rz(0.65027872) q[1];
sx q[1];
rz(-1.283409) q[1];
sx q[1];
rz(0.78279701) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37047526) q[0];
sx q[0];
rz(-1.820571) q[0];
sx q[0];
rz(-3.0779548) q[0];
x q[1];
rz(1.2230258) q[2];
sx q[2];
rz(-2.7070621) q[2];
sx q[2];
rz(-2.6968616) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70071917) q[1];
sx q[1];
rz(-2.6881725) q[1];
sx q[1];
rz(-1.8666301) q[1];
rz(-pi) q[2];
rz(0.31580117) q[3];
sx q[3];
rz(-1.6557906) q[3];
sx q[3];
rz(1.0335361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.89447442) q[2];
sx q[2];
rz(-1.0049745) q[2];
sx q[2];
rz(-3.0257814) q[2];
rz(1.5995021) q[3];
sx q[3];
rz(-0.096297979) q[3];
sx q[3];
rz(1.0533062) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88749921) q[0];
sx q[0];
rz(-2.5920581) q[0];
sx q[0];
rz(0.19533531) q[0];
rz(-0.37503606) q[1];
sx q[1];
rz(-1.6655567) q[1];
sx q[1];
rz(0.23981747) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6759684) q[0];
sx q[0];
rz(-1.3249319) q[0];
sx q[0];
rz(1.4093536) q[0];
rz(1.8717143) q[2];
sx q[2];
rz(-1.7728724) q[2];
sx q[2];
rz(-0.84504715) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6310196) q[1];
sx q[1];
rz(-2.3595855) q[1];
sx q[1];
rz(1.4960947) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2858743) q[3];
sx q[3];
rz(-3.120003) q[3];
sx q[3];
rz(1.3947226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.6372765) q[2];
sx q[2];
rz(-0.80233032) q[2];
sx q[2];
rz(1.3298539) q[2];
rz(1.3416393) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(-1.8168861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(0.66816107) q[0];
rz(1.6502624) q[1];
sx q[1];
rz(-0.69258339) q[1];
sx q[1];
rz(1.0659165) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50915584) q[0];
sx q[0];
rz(-1.8189438) q[0];
sx q[0];
rz(0.042225348) q[0];
x q[1];
rz(-2.9124444) q[2];
sx q[2];
rz(-0.56729588) q[2];
sx q[2];
rz(2.0958401) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.26820688) q[1];
sx q[1];
rz(-1.2899439) q[1];
sx q[1];
rz(0.98209776) q[1];
rz(-pi) q[2];
rz(-2.2320281) q[3];
sx q[3];
rz(-1.4759016) q[3];
sx q[3];
rz(-0.50300099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9540017) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(-3.1090453) q[2];
rz(2.7815946) q[3];
sx q[3];
rz(-2.0149752) q[3];
sx q[3];
rz(0.79157296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86768326) q[0];
sx q[0];
rz(-1.5554579) q[0];
sx q[0];
rz(-0.70670635) q[0];
rz(1.2061521) q[1];
sx q[1];
rz(-0.34148347) q[1];
sx q[1];
rz(1.6548086) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33071163) q[0];
sx q[0];
rz(-1.8510305) q[0];
sx q[0];
rz(2.7964554) q[0];
x q[1];
rz(0.47866486) q[2];
sx q[2];
rz(-2.475111) q[2];
sx q[2];
rz(-1.1043617) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.023991) q[1];
sx q[1];
rz(-1.7204493) q[1];
sx q[1];
rz(0.28848044) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3151822) q[3];
sx q[3];
rz(-1.5757757) q[3];
sx q[3];
rz(1.6299562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5965745) q[2];
sx q[2];
rz(-2.2480965) q[2];
sx q[2];
rz(2.13307) q[2];
rz(-2.0452943) q[3];
sx q[3];
rz(-1.9129646) q[3];
sx q[3];
rz(-0.1299468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0356045) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(-1.6812356) q[0];
rz(1.5530855) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(-0.016074093) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.21503) q[0];
sx q[0];
rz(-2.3030871) q[0];
sx q[0];
rz(0.82920427) q[0];
rz(0.0083382567) q[2];
sx q[2];
rz(-2.001488) q[2];
sx q[2];
rz(-0.24239937) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4453455) q[1];
sx q[1];
rz(-2.1228585) q[1];
sx q[1];
rz(-0.50260431) q[1];
rz(-1.5752951) q[3];
sx q[3];
rz(-1.2084949) q[3];
sx q[3];
rz(0.11210657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0323223) q[2];
sx q[2];
rz(-2.1123999) q[2];
sx q[2];
rz(2.5197022) q[2];
rz(2.0444929) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(2.9212852) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19875232) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(-2.2348256) q[0];
rz(-0.81470195) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(-1.9168568) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2621433) q[0];
sx q[0];
rz(-2.1855542) q[0];
sx q[0];
rz(0.26457796) q[0];
x q[1];
rz(-1.1512418) q[2];
sx q[2];
rz(-2.662979) q[2];
sx q[2];
rz(1.4065557) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3983706) q[1];
sx q[1];
rz(-2.1757158) q[1];
sx q[1];
rz(-0.066992316) q[1];
rz(-pi) q[2];
rz(2.740432) q[3];
sx q[3];
rz(-1.6725371) q[3];
sx q[3];
rz(1.7743558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1427052) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(-0.93969807) q[2];
rz(-0.21329221) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(-1.7512158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1619103) q[0];
sx q[0];
rz(-2.1770711) q[0];
sx q[0];
rz(-2.5612223) q[0];
rz(2.0866108) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(0.70077983) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38802559) q[0];
sx q[0];
rz(-1.1209079) q[0];
sx q[0];
rz(0.1247503) q[0];
rz(-1.6875661) q[2];
sx q[2];
rz(-1.2532557) q[2];
sx q[2];
rz(-1.5460528) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.82579457) q[1];
sx q[1];
rz(-1.0072395) q[1];
sx q[1];
rz(-3.0347996) q[1];
x q[2];
rz(-0.84007646) q[3];
sx q[3];
rz(-1.6335765) q[3];
sx q[3];
rz(-1.1508133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7769527) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(-0.022162612) q[2];
rz(0.74470216) q[3];
sx q[3];
rz(-2.0245168) q[3];
sx q[3];
rz(-0.40063342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78684029) q[0];
sx q[0];
rz(-2.1180034) q[0];
sx q[0];
rz(1.4165075) q[0];
rz(-1.7658866) q[1];
sx q[1];
rz(-1.7388652) q[1];
sx q[1];
rz(-2.2498806) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4448924) q[0];
sx q[0];
rz(-2.3377044) q[0];
sx q[0];
rz(-2.9219887) q[0];
rz(2.5052091) q[2];
sx q[2];
rz(-0.6837662) q[2];
sx q[2];
rz(0.46905876) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2878694) q[1];
sx q[1];
rz(-0.81765491) q[1];
sx q[1];
rz(2.5649648) q[1];
rz(-1.9705087) q[3];
sx q[3];
rz(-0.69972316) q[3];
sx q[3];
rz(-1.2801998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6531758) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(-0.78424224) q[2];
rz(2.6358321) q[3];
sx q[3];
rz(-0.85251802) q[3];
sx q[3];
rz(2.908356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6655675) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(2.419557) q[0];
rz(-0.33323914) q[1];
sx q[1];
rz(-1.9457341) q[1];
sx q[1];
rz(-1.3649712) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1418755) q[0];
sx q[0];
rz(-1.5293855) q[0];
sx q[0];
rz(2.3977604) q[0];
rz(0.083655595) q[2];
sx q[2];
rz(-1.9803279) q[2];
sx q[2];
rz(-1.6905897) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.078799876) q[1];
sx q[1];
rz(-1.2275057) q[1];
sx q[1];
rz(0.72453665) q[1];
x q[2];
rz(0.19946675) q[3];
sx q[3];
rz(-1.27928) q[3];
sx q[3];
rz(-0.22649543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0329131) q[2];
sx q[2];
rz(-1.379517) q[2];
sx q[2];
rz(1.2314679) q[2];
rz(0.03406295) q[3];
sx q[3];
rz(-1.27682) q[3];
sx q[3];
rz(-2.506822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(2.05474) q[0];
sx q[0];
rz(-0.56607902) q[0];
sx q[0];
rz(1.4779133) q[0];
rz(2.058303) q[1];
sx q[1];
rz(-1.7419086) q[1];
sx q[1];
rz(2.1733984) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4817081) q[0];
sx q[0];
rz(-1.9511001) q[0];
sx q[0];
rz(-0.36614059) q[0];
rz(-pi) q[1];
rz(0.88629006) q[2];
sx q[2];
rz(-0.61294014) q[2];
sx q[2];
rz(2.431589) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8755175) q[1];
sx q[1];
rz(-2.7873758) q[1];
sx q[1];
rz(-0.43550272) q[1];
rz(-pi) q[2];
rz(-1.5714985) q[3];
sx q[3];
rz(-0.32214221) q[3];
sx q[3];
rz(0.57442564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0633462) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(3.1402804) q[2];
rz(1.1408268) q[3];
sx q[3];
rz(-1.3169378) q[3];
sx q[3];
rz(-1.3425945) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7447727) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(-2.7753579) q[1];
sx q[1];
rz(-1.9402505) q[1];
sx q[1];
rz(-1.8016626) q[1];
rz(0.33404074) q[2];
sx q[2];
rz(-2.3636849) q[2];
sx q[2];
rz(-1.452527) q[2];
rz(1.0127388) q[3];
sx q[3];
rz(-1.8752718) q[3];
sx q[3];
rz(0.74753052) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];