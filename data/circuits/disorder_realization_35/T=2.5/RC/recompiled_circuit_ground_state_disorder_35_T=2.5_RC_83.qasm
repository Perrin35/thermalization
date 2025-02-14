OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2090981) q[0];
sx q[0];
rz(0.068109186) q[0];
sx q[0];
rz(13.343233) q[0];
rz(1.8344954) q[1];
sx q[1];
rz(-2.1721462) q[1];
sx q[1];
rz(-1.3219272) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711219) q[0];
sx q[0];
rz(-1.4323455) q[0];
sx q[0];
rz(-2.2884503) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49643699) q[2];
sx q[2];
rz(-1.8775216) q[2];
sx q[2];
rz(2.6076743) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2474587) q[1];
sx q[1];
rz(-1.9795291) q[1];
sx q[1];
rz(0.031886727) q[1];
x q[2];
rz(-0.98422756) q[3];
sx q[3];
rz(-2.7968458) q[3];
sx q[3];
rz(-2.1380827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0778568) q[2];
sx q[2];
rz(-1.7813762) q[2];
sx q[2];
rz(0.35749164) q[2];
rz(0.56719559) q[3];
sx q[3];
rz(-1.7287858) q[3];
sx q[3];
rz(0.43958694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95328632) q[0];
sx q[0];
rz(-2.8585241) q[0];
sx q[0];
rz(-1.8081007) q[0];
rz(-2.7045344) q[1];
sx q[1];
rz(-2.5407365) q[1];
sx q[1];
rz(0.83998799) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49522049) q[0];
sx q[0];
rz(-2.6073808) q[0];
sx q[0];
rz(-1.87508) q[0];
x q[1];
rz(-0.2510906) q[2];
sx q[2];
rz(-2.5320964) q[2];
sx q[2];
rz(1.9025546) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.67761356) q[1];
sx q[1];
rz(-0.85444728) q[1];
sx q[1];
rz(-0.78732995) q[1];
rz(-pi) q[2];
rz(1.6343104) q[3];
sx q[3];
rz(-2.5908785) q[3];
sx q[3];
rz(0.39001071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4240894) q[2];
sx q[2];
rz(-1.4428029) q[2];
sx q[2];
rz(1.7355512) q[2];
rz(0.16215912) q[3];
sx q[3];
rz(-1.5608965) q[3];
sx q[3];
rz(-0.46419188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9185987) q[0];
sx q[0];
rz(-3.0610237) q[0];
sx q[0];
rz(2.719847) q[0];
rz(2.2987507) q[1];
sx q[1];
rz(-1.755736) q[1];
sx q[1];
rz(-0.27580321) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11711794) q[0];
sx q[0];
rz(-1.3399276) q[0];
sx q[0];
rz(2.8083715) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7925198) q[2];
sx q[2];
rz(-1.1100551) q[2];
sx q[2];
rz(2.0784476) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0581857) q[1];
sx q[1];
rz(-2.8166933) q[1];
sx q[1];
rz(2.9778381) q[1];
x q[2];
rz(-3.0432947) q[3];
sx q[3];
rz(-2.2701828) q[3];
sx q[3];
rz(-0.15546945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.80321035) q[2];
sx q[2];
rz(-0.79757491) q[2];
sx q[2];
rz(-2.4513643) q[2];
rz(0.35987443) q[3];
sx q[3];
rz(-1.3709603) q[3];
sx q[3];
rz(-2.7580822) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91367078) q[0];
sx q[0];
rz(-2.0991195) q[0];
sx q[0];
rz(-0.87164718) q[0];
rz(-0.40920416) q[1];
sx q[1];
rz(-2.6458461) q[1];
sx q[1];
rz(-2.8292378) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46074122) q[0];
sx q[0];
rz(-0.4609403) q[0];
sx q[0];
rz(0.26682202) q[0];
rz(2.6957507) q[2];
sx q[2];
rz(-0.22805691) q[2];
sx q[2];
rz(-2.7331309) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.40939366) q[1];
sx q[1];
rz(-2.7229561) q[1];
sx q[1];
rz(2.6217209) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4931025) q[3];
sx q[3];
rz(-1.4146418) q[3];
sx q[3];
rz(0.79532571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.13545869) q[2];
sx q[2];
rz(-1.4486855) q[2];
sx q[2];
rz(2.3166166) q[2];
rz(-1.4247591) q[3];
sx q[3];
rz(-2.8202839) q[3];
sx q[3];
rz(2.7062374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.41998) q[0];
sx q[0];
rz(-1.5897607) q[0];
sx q[0];
rz(2.2671674) q[0];
rz(2.8127316) q[1];
sx q[1];
rz(-1.7560274) q[1];
sx q[1];
rz(-1.0985451) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4916353) q[0];
sx q[0];
rz(-1.1960317) q[0];
sx q[0];
rz(1.4521506) q[0];
rz(0.036089049) q[2];
sx q[2];
rz(-2.3054625) q[2];
sx q[2];
rz(0.058908894) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6381298) q[1];
sx q[1];
rz(-1.0833122) q[1];
sx q[1];
rz(2.8007228) q[1];
rz(-pi) q[2];
rz(0.27733488) q[3];
sx q[3];
rz(-0.99858741) q[3];
sx q[3];
rz(-2.8830143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.141779) q[2];
sx q[2];
rz(-1.9580611) q[2];
sx q[2];
rz(-0.67406526) q[2];
rz(-1.4874124) q[3];
sx q[3];
rz(-1.4451278) q[3];
sx q[3];
rz(2.8580247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0312626) q[0];
sx q[0];
rz(-2.1188348) q[0];
sx q[0];
rz(-1.7559825) q[0];
rz(-1.1194057) q[1];
sx q[1];
rz(-1.6740572) q[1];
sx q[1];
rz(1.3853692) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7840609) q[0];
sx q[0];
rz(-2.1936532) q[0];
sx q[0];
rz(-2.8297367) q[0];
x q[1];
rz(-2.643061) q[2];
sx q[2];
rz(-1.6462925) q[2];
sx q[2];
rz(0.99270644) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0138356) q[1];
sx q[1];
rz(-2.7633939) q[1];
sx q[1];
rz(-2.5950123) q[1];
rz(-1.5343666) q[3];
sx q[3];
rz(-1.4163989) q[3];
sx q[3];
rz(-1.2702219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0783656) q[2];
sx q[2];
rz(-2.207022) q[2];
sx q[2];
rz(0.37609491) q[2];
rz(-2.0070576) q[3];
sx q[3];
rz(-0.36342707) q[3];
sx q[3];
rz(-2.334972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019808708) q[0];
sx q[0];
rz(-0.60096318) q[0];
sx q[0];
rz(-0.10251775) q[0];
rz(0.58492297) q[1];
sx q[1];
rz(-0.94766098) q[1];
sx q[1];
rz(-1.1835416) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4393504) q[0];
sx q[0];
rz(-0.19225129) q[0];
sx q[0];
rz(0.73701145) q[0];
rz(-pi) q[1];
rz(2.7223515) q[2];
sx q[2];
rz(-0.70690522) q[2];
sx q[2];
rz(-2.8173994) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.592748) q[1];
sx q[1];
rz(-1.8251437) q[1];
sx q[1];
rz(0.91491048) q[1];
x q[2];
rz(0.46681301) q[3];
sx q[3];
rz(-1.3315505) q[3];
sx q[3];
rz(-1.5249263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3465603) q[2];
sx q[2];
rz(-2.2117895) q[2];
sx q[2];
rz(-0.19980508) q[2];
rz(2.0810769) q[3];
sx q[3];
rz(-2.6001055) q[3];
sx q[3];
rz(-3.0919891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.26315966) q[0];
sx q[0];
rz(-1.4819772) q[0];
sx q[0];
rz(-2.8286381) q[0];
rz(2.1921659) q[1];
sx q[1];
rz(-1.3525617) q[1];
sx q[1];
rz(-1.4535646) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4623186) q[0];
sx q[0];
rz(-1.3841713) q[0];
sx q[0];
rz(1.4786722) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5220479) q[2];
sx q[2];
rz(-1.6352374) q[2];
sx q[2];
rz(0.10945129) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.955757) q[1];
sx q[1];
rz(-1.794853) q[1];
sx q[1];
rz(-0.70835374) q[1];
x q[2];
rz(0.45519228) q[3];
sx q[3];
rz(-1.2641915) q[3];
sx q[3];
rz(2.2835177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.35385418) q[2];
sx q[2];
rz(-2.165803) q[2];
sx q[2];
rz(1.806951) q[2];
rz(1.3245964) q[3];
sx q[3];
rz(-0.66671222) q[3];
sx q[3];
rz(1.9273531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.1052833) q[0];
sx q[0];
rz(-2.5258625) q[0];
sx q[0];
rz(0.69806725) q[0];
rz(-2.2659194) q[1];
sx q[1];
rz(-1.061729) q[1];
sx q[1];
rz(2.7811513) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1669641) q[0];
sx q[0];
rz(-0.95050838) q[0];
sx q[0];
rz(-1.0247158) q[0];
rz(-2.7588604) q[2];
sx q[2];
rz(-1.6305411) q[2];
sx q[2];
rz(-2.003423) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1680582) q[1];
sx q[1];
rz(-2.413001) q[1];
sx q[1];
rz(2.161977) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8713851) q[3];
sx q[3];
rz(-2.5761309) q[3];
sx q[3];
rz(2.5266441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2299819) q[2];
sx q[2];
rz(-1.4260099) q[2];
sx q[2];
rz(-0.81725517) q[2];
rz(-0.83003712) q[3];
sx q[3];
rz(-0.54098141) q[3];
sx q[3];
rz(0.219492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.5068186) q[0];
sx q[0];
rz(-2.2932678) q[0];
sx q[0];
rz(-0.032055227) q[0];
rz(1.3196779) q[1];
sx q[1];
rz(-1.3751043) q[1];
sx q[1];
rz(-0.65418902) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35857707) q[0];
sx q[0];
rz(-0.23046432) q[0];
sx q[0];
rz(-0.5490659) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45126057) q[2];
sx q[2];
rz(-2.3006744) q[2];
sx q[2];
rz(-0.40715363) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4465877) q[1];
sx q[1];
rz(-0.69952974) q[1];
sx q[1];
rz(-2.3746731) q[1];
x q[2];
rz(-0.28367646) q[3];
sx q[3];
rz(-0.95796889) q[3];
sx q[3];
rz(-3.0949288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1303945) q[2];
sx q[2];
rz(-2.0989213) q[2];
sx q[2];
rz(1.0065669) q[2];
rz(2.448163) q[3];
sx q[3];
rz(-1.3207685) q[3];
sx q[3];
rz(-1.2815732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4257767) q[0];
sx q[0];
rz(-0.67831138) q[0];
sx q[0];
rz(-0.074445733) q[0];
rz(0.88059942) q[1];
sx q[1];
rz(-2.0279299) q[1];
sx q[1];
rz(0.50150064) q[1];
rz(0.51787372) q[2];
sx q[2];
rz(-0.93954115) q[2];
sx q[2];
rz(2.1672135) q[2];
rz(-0.69622688) q[3];
sx q[3];
rz(-1.8693557) q[3];
sx q[3];
rz(-2.1341677) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
