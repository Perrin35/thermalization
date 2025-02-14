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
rz(2.8528557) q[0];
sx q[0];
rz(-0.69536916) q[0];
sx q[0];
rz(-2.8737972) q[0];
rz(-2.7195622) q[1];
sx q[1];
rz(-0.92075092) q[1];
sx q[1];
rz(1.86778) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.081163) q[0];
sx q[0];
rz(-1.3295242) q[0];
sx q[0];
rz(-0.046242899) q[0];
rz(-pi) q[1];
rz(0.60611764) q[2];
sx q[2];
rz(-0.96783468) q[2];
sx q[2];
rz(-2.9204592) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5112524) q[1];
sx q[1];
rz(-0.66424886) q[1];
sx q[1];
rz(3.0306007) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7142322) q[3];
sx q[3];
rz(-1.535245) q[3];
sx q[3];
rz(1.2811023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5219118) q[2];
sx q[2];
rz(-0.64138594) q[2];
sx q[2];
rz(-0.042595159) q[2];
rz(-0.28111449) q[3];
sx q[3];
rz(-1.5646076) q[3];
sx q[3];
rz(2.5458096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5113145) q[0];
sx q[0];
rz(-1.9673286) q[0];
sx q[0];
rz(-1.0761155) q[0];
rz(2.1108744) q[1];
sx q[1];
rz(-1.1117671) q[1];
sx q[1];
rz(-1.8310742) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.977518) q[0];
sx q[0];
rz(-2.4830472) q[0];
sx q[0];
rz(-1.6048628) q[0];
rz(-pi) q[1];
rz(-0.4035455) q[2];
sx q[2];
rz(-0.87658823) q[2];
sx q[2];
rz(-2.780404) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0243822) q[1];
sx q[1];
rz(-1.4522933) q[1];
sx q[1];
rz(-2.973053) q[1];
x q[2];
rz(0.63046534) q[3];
sx q[3];
rz(-0.63172904) q[3];
sx q[3];
rz(0.28030685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.92404667) q[2];
sx q[2];
rz(-0.7520389) q[2];
sx q[2];
rz(2.1853866) q[2];
rz(1.8337967) q[3];
sx q[3];
rz(-2.3266413) q[3];
sx q[3];
rz(3.0202878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2995375) q[0];
sx q[0];
rz(-2.6326023) q[0];
sx q[0];
rz(-3.1053542) q[0];
rz(-0.80129519) q[1];
sx q[1];
rz(-1.5943297) q[1];
sx q[1];
rz(1.3005728) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7177082) q[0];
sx q[0];
rz(-2.9489655) q[0];
sx q[0];
rz(2.1602551) q[0];
rz(-pi) q[1];
rz(0.043134281) q[2];
sx q[2];
rz(-0.97402527) q[2];
sx q[2];
rz(-1.5297001) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.76884586) q[1];
sx q[1];
rz(-2.0018491) q[1];
sx q[1];
rz(-2.194904) q[1];
rz(-pi) q[2];
rz(1.9375422) q[3];
sx q[3];
rz(-0.91454187) q[3];
sx q[3];
rz(1.3964069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3761882) q[2];
sx q[2];
rz(-1.3724962) q[2];
sx q[2];
rz(-2.9236531) q[2];
rz(2.229522) q[3];
sx q[3];
rz(-1.6488766) q[3];
sx q[3];
rz(0.85165858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7243778) q[0];
sx q[0];
rz(-2.0700924) q[0];
sx q[0];
rz(-2.3260314) q[0];
rz(2.0888445) q[1];
sx q[1];
rz(-0.88200724) q[1];
sx q[1];
rz(1.7900593) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4690875) q[0];
sx q[0];
rz(-2.6658305) q[0];
sx q[0];
rz(-2.9722054) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6585245) q[2];
sx q[2];
rz(-0.85614877) q[2];
sx q[2];
rz(-1.8217979) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0114944) q[1];
sx q[1];
rz(-1.5421499) q[1];
sx q[1];
rz(-2.0462034) q[1];
rz(-2.0729321) q[3];
sx q[3];
rz(-0.91864139) q[3];
sx q[3];
rz(2.7606635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.012933) q[2];
sx q[2];
rz(-1.1933051) q[2];
sx q[2];
rz(-2.7093757) q[2];
rz(0.84248078) q[3];
sx q[3];
rz(-2.1079) q[3];
sx q[3];
rz(0.99561083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9812444) q[0];
sx q[0];
rz(-0.46537414) q[0];
sx q[0];
rz(0.55111849) q[0];
rz(2.5235858) q[1];
sx q[1];
rz(-1.2851241) q[1];
sx q[1];
rz(-1.0689703) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7525714) q[0];
sx q[0];
rz(-2.4271936) q[0];
sx q[0];
rz(2.8117489) q[0];
rz(-pi) q[1];
rz(-0.48241291) q[2];
sx q[2];
rz(-1.97792) q[2];
sx q[2];
rz(2.3633752) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3772419) q[1];
sx q[1];
rz(-0.83989401) q[1];
sx q[1];
rz(-2.9725171) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51286814) q[3];
sx q[3];
rz(-1.4445496) q[3];
sx q[3];
rz(0.083663705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3507639) q[2];
sx q[2];
rz(-0.5841693) q[2];
sx q[2];
rz(-2.186415) q[2];
rz(-2.5480934) q[3];
sx q[3];
rz(-2.1953526) q[3];
sx q[3];
rz(1.3957297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7668358) q[0];
sx q[0];
rz(-0.042348472) q[0];
sx q[0];
rz(-1.391885) q[0];
rz(-1.1426686) q[1];
sx q[1];
rz(-1.3418158) q[1];
sx q[1];
rz(-1.3124189) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6630733) q[0];
sx q[0];
rz(-2.0475302) q[0];
sx q[0];
rz(-0.8404503) q[0];
rz(-pi) q[1];
rz(0.67709558) q[2];
sx q[2];
rz(-0.99092084) q[2];
sx q[2];
rz(3.0486272) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7723263) q[1];
sx q[1];
rz(-1.6246987) q[1];
sx q[1];
rz(0.51477706) q[1];
x q[2];
rz(1.8862861) q[3];
sx q[3];
rz(-1.5157561) q[3];
sx q[3];
rz(2.6590527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.352508) q[2];
sx q[2];
rz(-2.3859873) q[2];
sx q[2];
rz(-1.4136723) q[2];
rz(2.6740354) q[3];
sx q[3];
rz(-0.65842015) q[3];
sx q[3];
rz(0.55268923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6308052) q[0];
sx q[0];
rz(-0.063022114) q[0];
sx q[0];
rz(-2.4600273) q[0];
rz(-1.956578) q[1];
sx q[1];
rz(-0.76528913) q[1];
sx q[1];
rz(2.3550745) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98726051) q[0];
sx q[0];
rz(-1.2367147) q[0];
sx q[0];
rz(0.41436899) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94633905) q[2];
sx q[2];
rz(-1.985744) q[2];
sx q[2];
rz(0.13042658) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29132641) q[1];
sx q[1];
rz(-0.89610142) q[1];
sx q[1];
rz(0.13723252) q[1];
rz(-1.5249355) q[3];
sx q[3];
rz(-2.4099382) q[3];
sx q[3];
rz(-1.589244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6098392) q[2];
sx q[2];
rz(-2.9015151) q[2];
sx q[2];
rz(-1.5717724) q[2];
rz(-1.2872559) q[3];
sx q[3];
rz(-1.6183034) q[3];
sx q[3];
rz(2.1985445) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.542273) q[0];
sx q[0];
rz(-2.4241408) q[0];
sx q[0];
rz(1.9532816) q[0];
rz(3.1207454) q[1];
sx q[1];
rz(-2.3884845) q[1];
sx q[1];
rz(2.5426224) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63510977) q[0];
sx q[0];
rz(-1.0671033) q[0];
sx q[0];
rz(-3.0290108) q[0];
rz(1.125419) q[2];
sx q[2];
rz(-2.3756304) q[2];
sx q[2];
rz(-0.19419032) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4074583) q[1];
sx q[1];
rz(-2.2770513) q[1];
sx q[1];
rz(2.4962673) q[1];
rz(-pi) q[2];
rz(2.3009572) q[3];
sx q[3];
rz(-1.908769) q[3];
sx q[3];
rz(-0.69608758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8591566) q[2];
sx q[2];
rz(-1.891529) q[2];
sx q[2];
rz(2.6263728) q[2];
rz(-0.53269261) q[3];
sx q[3];
rz(-1.0849413) q[3];
sx q[3];
rz(2.2044619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027997967) q[0];
sx q[0];
rz(-2.4585215) q[0];
sx q[0];
rz(0.37044507) q[0];
rz(-1.0796374) q[1];
sx q[1];
rz(-0.41079435) q[1];
sx q[1];
rz(-0.058578514) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81007018) q[0];
sx q[0];
rz(-2.5096623) q[0];
sx q[0];
rz(0.26946108) q[0];
rz(-pi) q[1];
rz(-1.3809105) q[2];
sx q[2];
rz(-0.23139628) q[2];
sx q[2];
rz(-2.2105202) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1778922) q[1];
sx q[1];
rz(-2.3517015) q[1];
sx q[1];
rz(-3.0067985) q[1];
rz(-pi) q[2];
rz(2.1982369) q[3];
sx q[3];
rz(-2.3010212) q[3];
sx q[3];
rz(2.3495903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2728682) q[2];
sx q[2];
rz(-0.803002) q[2];
sx q[2];
rz(2.8105984) q[2];
rz(0.013414772) q[3];
sx q[3];
rz(-2.1881073) q[3];
sx q[3];
rz(-2.0851871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57824221) q[0];
sx q[0];
rz(-0.42911068) q[0];
sx q[0];
rz(0.44813928) q[0];
rz(3.0478802) q[1];
sx q[1];
rz(-2.8401076) q[1];
sx q[1];
rz(1.2154481) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5543723) q[0];
sx q[0];
rz(-1.1379114) q[0];
sx q[0];
rz(-0.032399633) q[0];
rz(-pi) q[1];
rz(1.1555668) q[2];
sx q[2];
rz(-0.80249635) q[2];
sx q[2];
rz(1.323483) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3642973) q[1];
sx q[1];
rz(-2.5357995) q[1];
sx q[1];
rz(-2.0531027) q[1];
x q[2];
rz(1.1724505) q[3];
sx q[3];
rz(-1.6961548) q[3];
sx q[3];
rz(-2.0490058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.965968) q[2];
sx q[2];
rz(-1.5529239) q[2];
sx q[2];
rz(2.442181) q[2];
rz(-1.2753963) q[3];
sx q[3];
rz(-1.1558775) q[3];
sx q[3];
rz(-1.8705961) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70915837) q[0];
sx q[0];
rz(-2.2812738) q[0];
sx q[0];
rz(-1.9914837) q[0];
rz(1.746183) q[1];
sx q[1];
rz(-1.2184873) q[1];
sx q[1];
rz(1.7582735) q[1];
rz(-1.4906314) q[2];
sx q[2];
rz(-0.90654984) q[2];
sx q[2];
rz(-0.64150099) q[2];
rz(0.23308417) q[3];
sx q[3];
rz(-0.47745104) q[3];
sx q[3];
rz(0.61460809) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
