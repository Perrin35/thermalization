OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7528485) q[0];
sx q[0];
rz(-0.53628439) q[0];
sx q[0];
rz(-0.94776881) q[0];
rz(-1.3287969) q[1];
sx q[1];
rz(4.4089945) q[1];
sx q[1];
rz(10.452527) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42067895) q[0];
sx q[0];
rz(-1.5479711) q[0];
sx q[0];
rz(0.37585293) q[0];
x q[1];
rz(-2.920354) q[2];
sx q[2];
rz(-1.0916296) q[2];
sx q[2];
rz(-1.1269119) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7386908) q[1];
sx q[1];
rz(-1.1876904) q[1];
sx q[1];
rz(0.72523592) q[1];
x q[2];
rz(-2.0184228) q[3];
sx q[3];
rz(-2.1493704) q[3];
sx q[3];
rz(2.6920126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1296922) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(-1.0502846) q[2];
rz(-2.0283279) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(-0.068107001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0691836) q[0];
sx q[0];
rz(-1.2453112) q[0];
sx q[0];
rz(2.8438399) q[0];
rz(-2.521926) q[1];
sx q[1];
rz(-1.0071808) q[1];
sx q[1];
rz(1.108095) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1381123) q[0];
sx q[0];
rz(-1.3618042) q[0];
sx q[0];
rz(-0.94603993) q[0];
rz(-pi) q[1];
rz(-2.7116398) q[2];
sx q[2];
rz(-1.0064831) q[2];
sx q[2];
rz(-2.058409) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.121671) q[1];
sx q[1];
rz(-2.6056075) q[1];
sx q[1];
rz(0.74384816) q[1];
x q[2];
rz(-2.8262029) q[3];
sx q[3];
rz(-2.4139801) q[3];
sx q[3];
rz(-1.417158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6790598) q[2];
sx q[2];
rz(-0.97683895) q[2];
sx q[2];
rz(-2.1662946) q[2];
rz(0.9179999) q[3];
sx q[3];
rz(-1.8564329) q[3];
sx q[3];
rz(-0.29618922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
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
rz(-2.179203) q[0];
sx q[0];
rz(-2.2900892) q[0];
sx q[0];
rz(-0.54291022) q[0];
rz(-2.2593598) q[1];
sx q[1];
rz(-1.135332) q[1];
sx q[1];
rz(0.96484819) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99795139) q[0];
sx q[0];
rz(-0.43295857) q[0];
sx q[0];
rz(-0.71822449) q[0];
x q[1];
rz(2.8446571) q[2];
sx q[2];
rz(-1.529971) q[2];
sx q[2];
rz(-0.21619913) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.112903) q[1];
sx q[1];
rz(-1.6351055) q[1];
sx q[1];
rz(1.6921922) q[1];
x q[2];
rz(-1.996702) q[3];
sx q[3];
rz(-2.689038) q[3];
sx q[3];
rz(1.2061314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.09459153) q[2];
sx q[2];
rz(-2.5307405) q[2];
sx q[2];
rz(-2.0084521) q[2];
rz(0.23162332) q[3];
sx q[3];
rz(-1.2730205) q[3];
sx q[3];
rz(-0.75727063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27947458) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(1.3765155) q[0];
rz(-2.6230985) q[1];
sx q[1];
rz(-1.8771749) q[1];
sx q[1];
rz(-2.8994697) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0644558) q[0];
sx q[0];
rz(-0.69364871) q[0];
sx q[0];
rz(0.25607381) q[0];
rz(0.02853407) q[2];
sx q[2];
rz(-2.7118073) q[2];
sx q[2];
rz(-2.1536749) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.081269216) q[1];
sx q[1];
rz(-0.4821061) q[1];
sx q[1];
rz(0.37270765) q[1];
rz(1.4297156) q[3];
sx q[3];
rz(-1.1564848) q[3];
sx q[3];
rz(-2.8431818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.018192856) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(-0.094853178) q[2];
rz(1.799396) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(0.43911394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577394) q[0];
sx q[0];
rz(-1.8308324) q[0];
sx q[0];
rz(-2.7807996) q[0];
rz(-1.7533253) q[1];
sx q[1];
rz(-1.8107982) q[1];
sx q[1];
rz(-2.0070019) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9903588) q[0];
sx q[0];
rz(-0.84999527) q[0];
sx q[0];
rz(2.4609341) q[0];
x q[1];
rz(0.88419948) q[2];
sx q[2];
rz(-2.6066385) q[2];
sx q[2];
rz(1.3605489) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0134296) q[1];
sx q[1];
rz(-2.014782) q[1];
sx q[1];
rz(1.7432937) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4868823) q[3];
sx q[3];
rz(-2.0062371) q[3];
sx q[3];
rz(-2.8403789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1363042) q[2];
sx q[2];
rz(-1.7207928) q[2];
sx q[2];
rz(2.5689382) q[2];
rz(-2.2128361) q[3];
sx q[3];
rz(-2.6199665) q[3];
sx q[3];
rz(1.9740392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8144433) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(-1.2493398) q[0];
rz(1.2231187) q[1];
sx q[1];
rz(-1.616281) q[1];
sx q[1];
rz(1.9893601) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5778225) q[0];
sx q[0];
rz(-1.4655136) q[0];
sx q[0];
rz(0.012751243) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9828898) q[2];
sx q[2];
rz(-2.0372143) q[2];
sx q[2];
rz(-2.7127624) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.17105477) q[1];
sx q[1];
rz(-1.4661403) q[1];
sx q[1];
rz(2.9403951) q[1];
rz(0.96529393) q[3];
sx q[3];
rz(-2.1552342) q[3];
sx q[3];
rz(-0.30129978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46889177) q[2];
sx q[2];
rz(-1.9273309) q[2];
sx q[2];
rz(2.1506298) q[2];
rz(-0.64783603) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(0.34753862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.25794849) q[0];
sx q[0];
rz(-0.22709665) q[0];
sx q[0];
rz(0.062285475) q[0];
rz(0.1858055) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(0.3947765) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5236854) q[0];
sx q[0];
rz(-0.46645188) q[0];
sx q[0];
rz(-0.54752366) q[0];
rz(-pi) q[1];
rz(1.8918745) q[2];
sx q[2];
rz(-1.2121388) q[2];
sx q[2];
rz(-1.1527588) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9638622) q[1];
sx q[1];
rz(-1.9964255) q[1];
sx q[1];
rz(-1.8427986) q[1];
rz(-0.23468252) q[3];
sx q[3];
rz(-0.15771401) q[3];
sx q[3];
rz(1.371067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.64849598) q[2];
sx q[2];
rz(-0.33005565) q[2];
sx q[2];
rz(0.27302343) q[2];
rz(1.3027044) q[3];
sx q[3];
rz(-1.3132934) q[3];
sx q[3];
rz(-2.8295512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39772314) q[0];
sx q[0];
rz(-1.1345154) q[0];
sx q[0];
rz(1.8564818) q[0];
rz(-1.6400736) q[1];
sx q[1];
rz(-1.39095) q[1];
sx q[1];
rz(-1.8008908) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2368603) q[0];
sx q[0];
rz(-1.245984) q[0];
sx q[0];
rz(-1.1572641) q[0];
x q[1];
rz(2.4457744) q[2];
sx q[2];
rz(-1.5317481) q[2];
sx q[2];
rz(-2.2829636) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8074639) q[1];
sx q[1];
rz(-2.6263116) q[1];
sx q[1];
rz(0.79828243) q[1];
rz(-pi) q[2];
rz(-2.826346) q[3];
sx q[3];
rz(-1.9936221) q[3];
sx q[3];
rz(1.3413615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1061219) q[2];
sx q[2];
rz(-0.80646986) q[2];
sx q[2];
rz(2.1179874) q[2];
rz(2.9566531) q[3];
sx q[3];
rz(-2.7513294) q[3];
sx q[3];
rz(-0.24188724) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0446562) q[0];
sx q[0];
rz(-2.1450295) q[0];
sx q[0];
rz(-1.6212844) q[0];
rz(-0.3301436) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(-2.3044589) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9595813) q[0];
sx q[0];
rz(-0.89996979) q[0];
sx q[0];
rz(-1.3120033) q[0];
rz(-1.8673973) q[2];
sx q[2];
rz(-2.1269848) q[2];
sx q[2];
rz(1.1013168) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3763807) q[1];
sx q[1];
rz(-1.4188671) q[1];
sx q[1];
rz(3.0357749) q[1];
rz(-pi) q[2];
rz(-2.0666396) q[3];
sx q[3];
rz(-1.3131724) q[3];
sx q[3];
rz(0.50592929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.60124406) q[2];
sx q[2];
rz(-2.0566172) q[2];
sx q[2];
rz(-0.17962757) q[2];
rz(-0.99572292) q[3];
sx q[3];
rz(-1.8928173) q[3];
sx q[3];
rz(-1.3109591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3778465) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(2.0843704) q[0];
rz(-0.10748848) q[1];
sx q[1];
rz(-1.8881533) q[1];
sx q[1];
rz(0.9799788) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1996795) q[0];
sx q[0];
rz(-1.4451888) q[0];
sx q[0];
rz(3.0903387) q[0];
rz(0.51771848) q[2];
sx q[2];
rz(-1.1283518) q[2];
sx q[2];
rz(-0.25416086) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.53303888) q[1];
sx q[1];
rz(-1.8095008) q[1];
sx q[1];
rz(2.0134316) q[1];
rz(-1.8633217) q[3];
sx q[3];
rz(-2.8046126) q[3];
sx q[3];
rz(-1.302759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8367299) q[2];
sx q[2];
rz(-1.2048081) q[2];
sx q[2];
rz(1.0277964) q[2];
rz(1.7547539) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(-0.28361472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86823157) q[0];
sx q[0];
rz(-1.0366476) q[0];
sx q[0];
rz(-1.114053) q[0];
rz(-0.83203075) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(2.1297395) q[2];
sx q[2];
rz(-1.4197299) q[2];
sx q[2];
rz(1.1001669) q[2];
rz(0.63315331) q[3];
sx q[3];
rz(-0.59637759) q[3];
sx q[3];
rz(-0.37806088) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
