OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.40795657) q[0];
sx q[0];
rz(-2.9562558) q[0];
sx q[0];
rz(1.6428525) q[0];
rz(2.9474131) q[1];
sx q[1];
rz(-2.7979538) q[1];
sx q[1];
rz(-2.763881) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8728589) q[0];
sx q[0];
rz(-1.3493363) q[0];
sx q[0];
rz(0.24747699) q[0];
x q[1];
rz(2.7201817) q[2];
sx q[2];
rz(-1.974989) q[2];
sx q[2];
rz(0.39502783) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7491748) q[1];
sx q[1];
rz(-2.2274744) q[1];
sx q[1];
rz(1.5276457) q[1];
x q[2];
rz(-1.463428) q[3];
sx q[3];
rz(-0.26623785) q[3];
sx q[3];
rz(3.0913946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7517884) q[2];
sx q[2];
rz(-1.7916388) q[2];
sx q[2];
rz(-2.8765615) q[2];
rz(-2.7122279) q[3];
sx q[3];
rz(-1.8931484) q[3];
sx q[3];
rz(-0.00092367729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4877743) q[0];
sx q[0];
rz(-2.9980897) q[0];
sx q[0];
rz(1.300746) q[0];
rz(0.067151345) q[1];
sx q[1];
rz(-0.90922272) q[1];
sx q[1];
rz(-2.2815509) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22500998) q[0];
sx q[0];
rz(-2.2448178) q[0];
sx q[0];
rz(-2.3300578) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4323787) q[2];
sx q[2];
rz(-1.0150036) q[2];
sx q[2];
rz(-1.9947987) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7262267) q[1];
sx q[1];
rz(-2.4176855) q[1];
sx q[1];
rz(2.2103146) q[1];
rz(-2.9577291) q[3];
sx q[3];
rz(-1.2753092) q[3];
sx q[3];
rz(-1.9744524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5106875) q[2];
sx q[2];
rz(-0.61669934) q[2];
sx q[2];
rz(-2.5751233) q[2];
rz(-2.412292) q[3];
sx q[3];
rz(-0.84435487) q[3];
sx q[3];
rz(2.9384889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1931964) q[0];
sx q[0];
rz(-1.0862792) q[0];
sx q[0];
rz(2.7753944) q[0];
rz(-1.5858448) q[1];
sx q[1];
rz(-1.0428753) q[1];
sx q[1];
rz(-0.050447024) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6696861) q[0];
sx q[0];
rz(-1.6499053) q[0];
sx q[0];
rz(-0.0418274) q[0];
rz(-pi) q[1];
rz(2.453746) q[2];
sx q[2];
rz(-1.9645994) q[2];
sx q[2];
rz(1.1325768) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.37990268) q[1];
sx q[1];
rz(-2.8643048) q[1];
sx q[1];
rz(-2.4891341) q[1];
rz(-pi) q[2];
rz(-1.1459703) q[3];
sx q[3];
rz(-1.7033938) q[3];
sx q[3];
rz(-1.6679004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.067165) q[2];
sx q[2];
rz(-1.0600435) q[2];
sx q[2];
rz(-1.830843) q[2];
rz(-1.358076) q[3];
sx q[3];
rz(-2.570593) q[3];
sx q[3];
rz(1.3091492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8686789) q[0];
sx q[0];
rz(-2.9153115) q[0];
sx q[0];
rz(-1.7568461) q[0];
rz(-1.9028496) q[1];
sx q[1];
rz(-1.2307931) q[1];
sx q[1];
rz(1.1741656) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096701972) q[0];
sx q[0];
rz(-2.6706123) q[0];
sx q[0];
rz(-1.1440008) q[0];
x q[1];
rz(-2.410589) q[2];
sx q[2];
rz(-2.7766697) q[2];
sx q[2];
rz(-0.66563767) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8983072) q[1];
sx q[1];
rz(-1.8970058) q[1];
sx q[1];
rz(1.0350569) q[1];
rz(-pi) q[2];
rz(2.6888108) q[3];
sx q[3];
rz(-2.2119129) q[3];
sx q[3];
rz(-0.52911093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43526402) q[2];
sx q[2];
rz(-1.3202983) q[2];
sx q[2];
rz(3.1026133) q[2];
rz(0.6428166) q[3];
sx q[3];
rz(-1.0182074) q[3];
sx q[3];
rz(-2.5207998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83127999) q[0];
sx q[0];
rz(-2.1273002) q[0];
sx q[0];
rz(0.020462791) q[0];
rz(-0.19275716) q[1];
sx q[1];
rz(-2.5242476) q[1];
sx q[1];
rz(-1.1654759) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11949018) q[0];
sx q[0];
rz(-2.179232) q[0];
sx q[0];
rz(-1.0761989) q[0];
rz(-pi) q[1];
rz(-0.54372676) q[2];
sx q[2];
rz(-2.3909878) q[2];
sx q[2];
rz(3.1341396) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9679035) q[1];
sx q[1];
rz(-0.5040796) q[1];
sx q[1];
rz(1.5289115) q[1];
rz(1.9475157) q[3];
sx q[3];
rz(-2.93749) q[3];
sx q[3];
rz(0.78777504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.3141025) q[2];
sx q[2];
rz(-0.9919439) q[2];
sx q[2];
rz(-2.6143383) q[2];
rz(0.54316795) q[3];
sx q[3];
rz(-0.68734622) q[3];
sx q[3];
rz(0.34272042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0602033) q[0];
sx q[0];
rz(-0.15601604) q[0];
sx q[0];
rz(-2.7348837) q[0];
rz(1.1609062) q[1];
sx q[1];
rz(-2.2229767) q[1];
sx q[1];
rz(-1.0103753) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5766616) q[0];
sx q[0];
rz(-1.0258249) q[0];
sx q[0];
rz(-0.14133639) q[0];
rz(-0.40596227) q[2];
sx q[2];
rz(-2.7926707) q[2];
sx q[2];
rz(2.5418848) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0370673) q[1];
sx q[1];
rz(-0.43018331) q[1];
sx q[1];
rz(0.52806116) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25185092) q[3];
sx q[3];
rz(-0.82149071) q[3];
sx q[3];
rz(1.7918394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0668209) q[2];
sx q[2];
rz(-1.8596884) q[2];
sx q[2];
rz(-1.034896) q[2];
rz(-1.3567989) q[3];
sx q[3];
rz(-2.5467338) q[3];
sx q[3];
rz(-2.518173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.1739625) q[0];
sx q[0];
rz(-1.6169463) q[0];
sx q[0];
rz(-2.8136643) q[0];
rz(-3.0294042) q[1];
sx q[1];
rz(-1.9393549) q[1];
sx q[1];
rz(-0.97253886) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5060726) q[0];
sx q[0];
rz(-2.1284261) q[0];
sx q[0];
rz(-2.7235759) q[0];
rz(2.8424524) q[2];
sx q[2];
rz(-1.1465596) q[2];
sx q[2];
rz(-1.5787293) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.211208) q[1];
sx q[1];
rz(-2.5519785) q[1];
sx q[1];
rz(-3.1077887) q[1];
rz(-pi) q[2];
rz(-2.6554606) q[3];
sx q[3];
rz(-2.4555169) q[3];
sx q[3];
rz(3.0239575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.605725) q[2];
sx q[2];
rz(-0.87656993) q[2];
sx q[2];
rz(-2.4884339) q[2];
rz(-0.43290916) q[3];
sx q[3];
rz(-0.50986367) q[3];
sx q[3];
rz(-2.9292817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(1.9369649) q[0];
sx q[0];
rz(-0.84119868) q[0];
sx q[0];
rz(-0.2247819) q[0];
rz(-0.79554355) q[1];
sx q[1];
rz(-1.8276151) q[1];
sx q[1];
rz(2.0733817) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5607308) q[0];
sx q[0];
rz(-2.1819072) q[0];
sx q[0];
rz(-0.9471425) q[0];
rz(-2.3722367) q[2];
sx q[2];
rz(-2.0343668) q[2];
sx q[2];
rz(3.0887512) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.33419106) q[1];
sx q[1];
rz(-2.4535123) q[1];
sx q[1];
rz(-0.040236878) q[1];
x q[2];
rz(1.3490476) q[3];
sx q[3];
rz(-1.6283247) q[3];
sx q[3];
rz(-1.8309995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3855359) q[2];
sx q[2];
rz(-1.3827366) q[2];
sx q[2];
rz(0.63560152) q[2];
rz(-2.8086737) q[3];
sx q[3];
rz(-2.1300485) q[3];
sx q[3];
rz(-0.46271589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3391089) q[0];
sx q[0];
rz(-0.026263069) q[0];
sx q[0];
rz(0.12301692) q[0];
rz(1.3975551) q[1];
sx q[1];
rz(-1.3879644) q[1];
sx q[1];
rz(-2.3618598) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5194979) q[0];
sx q[0];
rz(-0.72185282) q[0];
sx q[0];
rz(1.9869542) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4559007) q[2];
sx q[2];
rz(-1.4567516) q[2];
sx q[2];
rz(-1.4742849) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.349682) q[1];
sx q[1];
rz(-1.9432707) q[1];
sx q[1];
rz(-2.9279686) q[1];
x q[2];
rz(1.6980096) q[3];
sx q[3];
rz(-1.6774584) q[3];
sx q[3];
rz(0.17388177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.66703779) q[2];
sx q[2];
rz(-1.8755308) q[2];
sx q[2];
rz(-3.0511268) q[2];
rz(0.41680923) q[3];
sx q[3];
rz(-0.79206812) q[3];
sx q[3];
rz(0.99378234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6592634) q[0];
sx q[0];
rz(-1.8195131) q[0];
sx q[0];
rz(3.0521159) q[0];
rz(0.80884519) q[1];
sx q[1];
rz(-0.50061148) q[1];
sx q[1];
rz(2.083875) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63454483) q[0];
sx q[0];
rz(-1.9814361) q[0];
sx q[0];
rz(-1.157874) q[0];
rz(-pi) q[1];
rz(1.5799149) q[2];
sx q[2];
rz(-2.5895725) q[2];
sx q[2];
rz(0.43285757) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32409975) q[1];
sx q[1];
rz(-1.9655394) q[1];
sx q[1];
rz(-2.7604719) q[1];
rz(-1.5461033) q[3];
sx q[3];
rz(-2.3273483) q[3];
sx q[3];
rz(0.3972185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7384501) q[2];
sx q[2];
rz(-0.53692997) q[2];
sx q[2];
rz(2.7682847) q[2];
rz(0.25660723) q[3];
sx q[3];
rz(-1.720287) q[3];
sx q[3];
rz(0.074450113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8858717) q[0];
sx q[0];
rz(-1.5789565) q[0];
sx q[0];
rz(-1.4174905) q[0];
rz(-1.7120842) q[1];
sx q[1];
rz(-0.34965873) q[1];
sx q[1];
rz(0.8082334) q[1];
rz(-1.5612372) q[2];
sx q[2];
rz(-2.0295967) q[2];
sx q[2];
rz(1.5461736) q[2];
rz(-3.1028845) q[3];
sx q[3];
rz(-2.1716272) q[3];
sx q[3];
rz(1.4089331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
