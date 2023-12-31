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
rz(2.6842527) q[1];
sx q[1];
rz(-2.1973124) q[1];
sx q[1];
rz(1.9231208) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96555644) q[0];
sx q[0];
rz(-2.3210225) q[0];
sx q[0];
rz(-1.5216989) q[0];
rz(-3.0604612) q[2];
sx q[2];
rz(-0.46617026) q[2];
sx q[2];
rz(-0.98035882) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4301181) q[1];
sx q[1];
rz(-2.4415486) q[1];
sx q[1];
rz(-2.3858566) q[1];
x q[2];
rz(-2.9524607) q[3];
sx q[3];
rz(-1.7889708) q[3];
sx q[3];
rz(1.8060341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5228287) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(0.5775601) q[2];
rz(1.1497568) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(-2.4770588) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1537271) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(-2.7541449) q[0];
rz(2.2024343) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(1.739025) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57099045) q[0];
sx q[0];
rz(-1.6002858) q[0];
sx q[0];
rz(1.5667314) q[0];
rz(-pi) q[1];
rz(-0.55949975) q[2];
sx q[2];
rz(-1.7569949) q[2];
sx q[2];
rz(-0.49636832) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2991997) q[1];
sx q[1];
rz(-1.0789011) q[1];
sx q[1];
rz(-2.030034) q[1];
rz(1.6658695) q[3];
sx q[3];
rz(-2.2007696) q[3];
sx q[3];
rz(0.41439393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7188321) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(0.31769162) q[2];
rz(2.9348532) q[3];
sx q[3];
rz(-0.59967774) q[3];
sx q[3];
rz(-0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7725672) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(1.7279708) q[0];
rz(-0.47779045) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(0.40107045) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0553592) q[0];
sx q[0];
rz(-1.2463352) q[0];
sx q[0];
rz(3.0326891) q[0];
x q[1];
rz(2.4972649) q[2];
sx q[2];
rz(-1.7228848) q[2];
sx q[2];
rz(-0.49567859) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8071825) q[1];
sx q[1];
rz(-1.5771598) q[1];
sx q[1];
rz(-0.73015726) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19604711) q[3];
sx q[3];
rz(-1.184433) q[3];
sx q[3];
rz(2.8783609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5473189) q[2];
sx q[2];
rz(-1.6197562) q[2];
sx q[2];
rz(0.55580124) q[2];
rz(-0.9764955) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34898409) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(1.4439616) q[0];
rz(1.6216888) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(0.25340432) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93638203) q[0];
sx q[0];
rz(-1.0316327) q[0];
sx q[0];
rz(0.11986952) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0850111) q[2];
sx q[2];
rz(-1.5554785) q[2];
sx q[2];
rz(2.3682396) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1052103) q[1];
sx q[1];
rz(-1.6670334) q[1];
sx q[1];
rz(-0.63016816) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6494135) q[3];
sx q[3];
rz(-1.489349) q[3];
sx q[3];
rz(-0.34929517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0306586) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(2.3542662) q[2];
rz(-0.91286719) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(-2.1319938) q[3];
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
rz(2.426429) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(-3.0786247) q[0];
rz(-3.0175623) q[1];
sx q[1];
rz(-2.3359559) q[1];
sx q[1];
rz(-2.6834992) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9034018) q[0];
sx q[0];
rz(-1.5731249) q[0];
sx q[0];
rz(2.6989614) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45256726) q[2];
sx q[2];
rz(-2.1308225) q[2];
sx q[2];
rz(-0.82816154) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8372247) q[1];
sx q[1];
rz(-0.75290426) q[1];
sx q[1];
rz(1.3992845) q[1];
rz(-pi) q[2];
rz(1.0701418) q[3];
sx q[3];
rz(-1.6541012) q[3];
sx q[3];
rz(-1.3798151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1725585) q[2];
sx q[2];
rz(-2.2183552) q[2];
sx q[2];
rz(0.22949533) q[2];
rz(3.138792) q[3];
sx q[3];
rz(-0.87001785) q[3];
sx q[3];
rz(-1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5620419) q[0];
sx q[0];
rz(-0.27844772) q[0];
sx q[0];
rz(-0.91947412) q[0];
rz(0.062285034) q[1];
sx q[1];
rz(-2.1376164) q[1];
sx q[1];
rz(1.2671635) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6309109) q[0];
sx q[0];
rz(-1.9598538) q[0];
sx q[0];
rz(-2.3657777) q[0];
rz(1.7153347) q[2];
sx q[2];
rz(-0.83206165) q[2];
sx q[2];
rz(1.9369672) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2515182) q[1];
sx q[1];
rz(-1.2444082) q[1];
sx q[1];
rz(2.9337487) q[1];
x q[2];
rz(-2.0165765) q[3];
sx q[3];
rz(-0.52281724) q[3];
sx q[3];
rz(1.0558053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.9617812) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(-2.5578257) q[2];
rz(-2.4328655) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(-0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07847438) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(1.0725347) q[0];
rz(-2.5947) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(-2.0297208) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5892964) q[0];
sx q[0];
rz(-2.6483905) q[0];
sx q[0];
rz(-1.769355) q[0];
rz(-pi) q[1];
rz(1.215559) q[2];
sx q[2];
rz(-1.0676427) q[2];
sx q[2];
rz(0.1728729) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.062473) q[1];
sx q[1];
rz(-1.8603431) q[1];
sx q[1];
rz(1.6505961) q[1];
rz(-pi) q[2];
rz(1.0681549) q[3];
sx q[3];
rz(-1.916269) q[3];
sx q[3];
rz(-0.84514602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.303858) q[2];
sx q[2];
rz(-1.6656817) q[2];
sx q[2];
rz(1.1676577) q[2];
rz(-1.6052823) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(-1.4060098) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4090356) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(2.4107966) q[0];
rz(-0.90019512) q[1];
sx q[1];
rz(-2.3370445) q[1];
sx q[1];
rz(-2.3866167) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6801493) q[0];
sx q[0];
rz(-1.7168683) q[0];
sx q[0];
rz(0.84772528) q[0];
rz(-pi) q[1];
x q[1];
rz(2.412699) q[2];
sx q[2];
rz(-2.0260099) q[2];
sx q[2];
rz(-0.58065562) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.054143993) q[1];
sx q[1];
rz(-0.99522299) q[1];
sx q[1];
rz(1.4364442) q[1];
rz(0.8927535) q[3];
sx q[3];
rz(-2.4340981) q[3];
sx q[3];
rz(1.0196109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.76453152) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(0.63684741) q[2];
rz(2.8751255) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(-1.5554957) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72717845) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(-1.138858) q[0];
rz(-0.75421929) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(-3.1220904) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4365387) q[0];
sx q[0];
rz(-2.9199613) q[0];
sx q[0];
rz(1.9428695) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9901864) q[2];
sx q[2];
rz(-1.3726241) q[2];
sx q[2];
rz(0.99934794) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.94432482) q[1];
sx q[1];
rz(-1.7515344) q[1];
sx q[1];
rz(-0.9428057) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9762514) q[3];
sx q[3];
rz(-1.8792361) q[3];
sx q[3];
rz(0.43499085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.0043682178) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(2.2926889) q[2];
rz(-2.7539339) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.6000115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7983109) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(-2.6570901) q[0];
rz(1.3867406) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(1.1482931) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023037993) q[0];
sx q[0];
rz(-1.430129) q[0];
sx q[0];
rz(1.592357) q[0];
rz(-pi) q[1];
rz(1.7933153) q[2];
sx q[2];
rz(-0.95138022) q[2];
sx q[2];
rz(-0.36048181) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82603589) q[1];
sx q[1];
rz(-0.6912187) q[1];
sx q[1];
rz(1.8408937) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6230574) q[3];
sx q[3];
rz(-1.2219547) q[3];
sx q[3];
rz(-2.7319752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2293573) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(-2.7764376) q[2];
rz(0.12864104) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(-0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1098332) q[0];
sx q[0];
rz(-0.84072996) q[0];
sx q[0];
rz(1.6054556) q[0];
rz(2.1784492) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(0.49114901) q[2];
sx q[2];
rz(-2.4158203) q[2];
sx q[2];
rz(-0.62266785) q[2];
rz(-2.1309489) q[3];
sx q[3];
rz(-1.5394566) q[3];
sx q[3];
rz(1.4395366) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
