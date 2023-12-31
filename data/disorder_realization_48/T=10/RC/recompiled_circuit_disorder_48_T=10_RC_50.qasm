OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.36800185) q[0];
sx q[0];
rz(-0.79080963) q[0];
sx q[0];
rz(-2.8074582) q[0];
rz(-0.45733991) q[1];
sx q[1];
rz(-0.9442803) q[1];
sx q[1];
rz(1.2184719) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63872913) q[0];
sx q[0];
rz(-1.5348866) q[0];
sx q[0];
rz(2.3907651) q[0];
rz(-pi) q[1];
rz(3.0604612) q[2];
sx q[2];
rz(-0.46617026) q[2];
sx q[2];
rz(0.98035882) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6005046) q[1];
sx q[1];
rz(-2.0588015) q[1];
sx q[1];
rz(-2.0946676) q[1];
rz(-pi) q[2];
rz(-1.7928042) q[3];
sx q[3];
rz(-1.7553925) q[3];
sx q[3];
rz(-2.8649462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5228287) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(0.5775601) q[2];
rz(-1.9918359) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.1537271) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(0.38744774) q[0];
rz(-0.93915835) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(-1.4025677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1419066) q[0];
sx q[0];
rz(-1.5748595) q[0];
sx q[0];
rz(0.02948972) q[0];
rz(-pi) q[1];
rz(-2.8005373) q[2];
sx q[2];
rz(-2.555072) q[2];
sx q[2];
rz(-1.3618493) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2991997) q[1];
sx q[1];
rz(-1.0789011) q[1];
sx q[1];
rz(2.030034) q[1];
rz(-pi) q[2];
rz(-2.5094633) q[3];
sx q[3];
rz(-1.4940133) q[3];
sx q[3];
rz(2.0413105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42276057) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(-2.823901) q[2];
rz(-0.20673949) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(1.3690255) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(1.4136219) q[0];
rz(-2.6638022) q[1];
sx q[1];
rz(-1.3505892) q[1];
sx q[1];
rz(0.40107045) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51940489) q[0];
sx q[0];
rz(-1.4675958) q[0];
sx q[0];
rz(1.2445356) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3813854) q[2];
sx q[2];
rz(-2.2064798) q[2];
sx q[2];
rz(0.96178255) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.89951) q[1];
sx q[1];
rz(-2.3009355) q[1];
sx q[1];
rz(1.5793369) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1242261) q[3];
sx q[3];
rz(-0.43101573) q[3];
sx q[3];
rz(-2.9197846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5473189) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(2.5857914) q[2];
rz(2.1650971) q[3];
sx q[3];
rz(-0.54978168) q[3];
sx q[3];
rz(-2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7926086) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(1.697631) q[0];
rz(-1.5199039) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(-2.8881883) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2052106) q[0];
sx q[0];
rz(-1.0316327) q[0];
sx q[0];
rz(3.0217231) q[0];
x q[1];
rz(-1.6035945) q[2];
sx q[2];
rz(-2.6555853) q[2];
sx q[2];
rz(-2.3151468) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.036382347) q[1];
sx q[1];
rz(-1.4745592) q[1];
sx q[1];
rz(-2.5114245) q[1];
rz(-pi) q[2];
rz(3.0074189) q[3];
sx q[3];
rz(-0.65376702) q[3];
sx q[3];
rz(1.8133481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.110934) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(2.3542662) q[2];
rz(0.91286719) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(-1.0095989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-2.426429) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(-0.062967904) q[0];
rz(-3.0175623) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(-0.45809349) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8040745) q[0];
sx q[0];
rz(-2.6989557) q[0];
sx q[0];
rz(-0.0054365693) q[0];
rz(-pi) q[1];
rz(0.96197084) q[2];
sx q[2];
rz(-1.95032) q[2];
sx q[2];
rz(-0.48987197) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3921567) q[1];
sx q[1];
rz(-1.6877618) q[1];
sx q[1];
rz(-2.3163296) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7430274) q[3];
sx q[3];
rz(-2.6346364) q[3];
sx q[3];
rz(2.7996922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9690341) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(0.22949533) q[2];
rz(0.0028006639) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(1.3389448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5620419) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(0.91947412) q[0];
rz(0.062285034) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(1.8744291) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7125268) q[0];
sx q[0];
rz(-2.292284) q[0];
sx q[0];
rz(-2.6119786) q[0];
rz(-0.15684776) q[2];
sx q[2];
rz(-2.3914797) q[2];
sx q[2];
rz(1.4175121) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.25176469) q[1];
sx q[1];
rz(-1.3740731) q[1];
sx q[1];
rz(-1.2377435) q[1];
rz(-0.24354981) q[3];
sx q[3];
rz(-2.0381513) q[3];
sx q[3];
rz(-1.5817643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.9617812) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(-0.58376694) q[2];
rz(-2.4328655) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(-0.023199737) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07847438) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(2.069058) q[0];
rz(2.5947) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(2.0297208) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5892964) q[0];
sx q[0];
rz(-2.6483905) q[0];
sx q[0];
rz(-1.769355) q[0];
x q[1];
rz(-0.53084897) q[2];
sx q[2];
rz(-1.8804272) q[2];
sx q[2];
rz(1.9206778) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.48549451) q[1];
sx q[1];
rz(-1.494325) q[1];
sx q[1];
rz(-2.8511726) q[1];
rz(-pi) q[2];
rz(-0.38970077) q[3];
sx q[3];
rz(-1.1003564) q[3];
sx q[3];
rz(0.54159347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.303858) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(1.973935) q[2];
rz(1.5363103) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(1.7355828) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7325571) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(2.4107966) q[0];
rz(0.90019512) q[1];
sx q[1];
rz(-2.3370445) q[1];
sx q[1];
rz(2.3866167) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6801493) q[0];
sx q[0];
rz(-1.7168683) q[0];
sx q[0];
rz(2.2938674) q[0];
rz(2.5078012) q[2];
sx q[2];
rz(-2.3049424) q[2];
sx q[2];
rz(-1.6939236) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6983812) q[1];
sx q[1];
rz(-1.6834007) q[1];
sx q[1];
rz(2.5618782) q[1];
rz(2.2488392) q[3];
sx q[3];
rz(-0.70749456) q[3];
sx q[3];
rz(1.0196109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3770611) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(2.5047452) q[2];
rz(2.8751255) q[3];
sx q[3];
rz(-1.0576495) q[3];
sx q[3];
rz(-1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72717845) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(1.138858) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(0.019502217) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4365387) q[0];
sx q[0];
rz(-0.22163135) q[0];
sx q[0];
rz(-1.1987232) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2150061) q[2];
sx q[2];
rz(-0.24878657) q[2];
sx q[2];
rz(-0.34005806) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1972678) q[1];
sx q[1];
rz(-1.7515344) q[1];
sx q[1];
rz(0.9428057) q[1];
rz(-2.0476258) q[3];
sx q[3];
rz(-0.34871021) q[3];
sx q[3];
rz(2.204012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.0043682178) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(0.84890378) q[2];
rz(-0.38765872) q[3];
sx q[3];
rz(-1.1281697) q[3];
sx q[3];
rz(-1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7983109) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(2.6570901) q[0];
rz(1.7548521) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(-1.1482931) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5507817) q[0];
sx q[0];
rz(-1.592144) q[0];
sx q[0];
rz(-3.000893) q[0];
rz(-pi) q[1];
rz(-2.5103288) q[2];
sx q[2];
rz(-1.3901276) q[2];
sx q[2];
rz(-1.0797015) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1710098) q[1];
sx q[1];
rz(-2.2323771) q[1];
sx q[1];
rz(0.2172825) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34928068) q[3];
sx q[3];
rz(-1.619907) q[3];
sx q[3];
rz(-1.1790566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2293573) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(-2.7764376) q[2];
rz(0.12864104) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(2.685759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0317595) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(-2.1784492) q[1];
sx q[1];
rz(-1.8704725) q[1];
sx q[1];
rz(2.0830547) q[1];
rz(1.1744432) q[2];
sx q[2];
rz(-2.195993) q[2];
sx q[2];
rz(3.1396951) q[2];
rz(1.5118602) q[3];
sx q[3];
rz(-0.56093506) q[3];
sx q[3];
rz(-0.18118071) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
