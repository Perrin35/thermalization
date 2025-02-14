OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.423288) q[0];
sx q[0];
rz(-2.365132) q[0];
sx q[0];
rz(-2.7756696) q[0];
rz(-2.794682) q[1];
sx q[1];
rz(3.4263098) q[1];
sx q[1];
rz(10.813536) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91639389) q[0];
sx q[0];
rz(-1.7728191) q[0];
sx q[0];
rz(1.5027769) q[0];
x q[1];
rz(0.69268815) q[2];
sx q[2];
rz(-1.3387623) q[2];
sx q[2];
rz(-1.2218066) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5947111) q[1];
sx q[1];
rz(-1.9597907) q[1];
sx q[1];
rz(1.0884029) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0991058) q[3];
sx q[3];
rz(-1.0431223) q[3];
sx q[3];
rz(1.8535525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1846788) q[2];
sx q[2];
rz(-0.92014402) q[2];
sx q[2];
rz(0.93519768) q[2];
rz(0.30185559) q[3];
sx q[3];
rz(-1.5993092) q[3];
sx q[3];
rz(2.7479318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.1644156) q[0];
sx q[0];
rz(-1.2657607) q[0];
sx q[0];
rz(0.49829495) q[0];
rz(1.7138819) q[1];
sx q[1];
rz(-1.2063113) q[1];
sx q[1];
rz(-2.0548342) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55288494) q[0];
sx q[0];
rz(-2.4037139) q[0];
sx q[0];
rz(2.7859475) q[0];
rz(1.837376) q[2];
sx q[2];
rz(-1.3836897) q[2];
sx q[2];
rz(-0.10548909) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9683428) q[1];
sx q[1];
rz(-1.1361215) q[1];
sx q[1];
rz(1.3938851) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8233775) q[3];
sx q[3];
rz(-1.5822142) q[3];
sx q[3];
rz(-2.5493357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.19865092) q[2];
sx q[2];
rz(-2.8309839) q[2];
sx q[2];
rz(-1.8093713) q[2];
rz(1.1474991) q[3];
sx q[3];
rz(-1.5787326) q[3];
sx q[3];
rz(-1.3326299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9320817) q[0];
sx q[0];
rz(-1.4028343) q[0];
sx q[0];
rz(-0.15723666) q[0];
rz(-1.9137742) q[1];
sx q[1];
rz(-0.20918748) q[1];
sx q[1];
rz(0.42207178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8764719) q[0];
sx q[0];
rz(-0.19959627) q[0];
sx q[0];
rz(1.100698) q[0];
rz(-1.3968758) q[2];
sx q[2];
rz(-1.2971767) q[2];
sx q[2];
rz(-1.5652986) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6032519) q[1];
sx q[1];
rz(-1.9515688) q[1];
sx q[1];
rz(3.046585) q[1];
rz(-2.3396569) q[3];
sx q[3];
rz(-1.7953331) q[3];
sx q[3];
rz(-0.79659407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.32855365) q[2];
sx q[2];
rz(-2.1812794) q[2];
sx q[2];
rz(-0.73406827) q[2];
rz(1.0775393) q[3];
sx q[3];
rz(-0.68818337) q[3];
sx q[3];
rz(0.075210007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5828736) q[0];
sx q[0];
rz(-1.0105157) q[0];
sx q[0];
rz(-0.016481312) q[0];
rz(0.074542848) q[1];
sx q[1];
rz(-2.6838979) q[1];
sx q[1];
rz(2.8864536) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1265565) q[0];
sx q[0];
rz(-1.714596) q[0];
sx q[0];
rz(2.6565927) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8066977) q[2];
sx q[2];
rz(-2.1112295) q[2];
sx q[2];
rz(1.188736) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0485766) q[1];
sx q[1];
rz(-2.3469791) q[1];
sx q[1];
rz(2.5203185) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8750853) q[3];
sx q[3];
rz(-0.36484584) q[3];
sx q[3];
rz(1.5675275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5269346) q[2];
sx q[2];
rz(-2.4444828) q[2];
sx q[2];
rz(-2.441791) q[2];
rz(0.091863306) q[3];
sx q[3];
rz(-1.9242761) q[3];
sx q[3];
rz(2.6868668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62234539) q[0];
sx q[0];
rz(-0.82038251) q[0];
sx q[0];
rz(-2.0297594) q[0];
rz(0.19019292) q[1];
sx q[1];
rz(-2.2564502) q[1];
sx q[1];
rz(-1.1598738) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1721508) q[0];
sx q[0];
rz(-1.5363201) q[0];
sx q[0];
rz(-1.9146754) q[0];
rz(-pi) q[1];
rz(2.9952107) q[2];
sx q[2];
rz(-2.2040151) q[2];
sx q[2];
rz(-0.13356471) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24007356) q[1];
sx q[1];
rz(-1.0283955) q[1];
sx q[1];
rz(1.2804922) q[1];
rz(-pi) q[2];
rz(-1.0489063) q[3];
sx q[3];
rz(-1.6097704) q[3];
sx q[3];
rz(-0.029595395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5708892) q[2];
sx q[2];
rz(-2.3184226) q[2];
sx q[2];
rz(2.9111351) q[2];
rz(-1.0620091) q[3];
sx q[3];
rz(-1.8657203) q[3];
sx q[3];
rz(-1.6331858) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45288169) q[0];
sx q[0];
rz(-1.0588366) q[0];
sx q[0];
rz(-1.7857312) q[0];
rz(0.52976766) q[1];
sx q[1];
rz(-2.3593088) q[1];
sx q[1];
rz(-1.9897602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065581948) q[0];
sx q[0];
rz(-1.9261203) q[0];
sx q[0];
rz(-1.4655345) q[0];
rz(-1.0291589) q[2];
sx q[2];
rz(-1.6525606) q[2];
sx q[2];
rz(1.077291) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.85688389) q[1];
sx q[1];
rz(-1.207106) q[1];
sx q[1];
rz(1.3631352) q[1];
x q[2];
rz(1.4532667) q[3];
sx q[3];
rz(-1.7835878) q[3];
sx q[3];
rz(0.26373395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6666224) q[2];
sx q[2];
rz(-0.64133659) q[2];
sx q[2];
rz(-0.46710157) q[2];
rz(-2.1304255) q[3];
sx q[3];
rz(-2.7979388) q[3];
sx q[3];
rz(-1.3694192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8491299) q[0];
sx q[0];
rz(-2.9865773) q[0];
sx q[0];
rz(-2.9169061) q[0];
rz(-0.16381964) q[1];
sx q[1];
rz(-2.8024709) q[1];
sx q[1];
rz(2.4952369) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40898973) q[0];
sx q[0];
rz(-1.8206795) q[0];
sx q[0];
rz(2.1134209) q[0];
x q[1];
rz(0.98585821) q[2];
sx q[2];
rz(-0.38160593) q[2];
sx q[2];
rz(1.3473036) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9440981) q[1];
sx q[1];
rz(-2.1131599) q[1];
sx q[1];
rz(1.4828862) q[1];
x q[2];
rz(-2.8180653) q[3];
sx q[3];
rz(-0.40641847) q[3];
sx q[3];
rz(-0.42662963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9859163) q[2];
sx q[2];
rz(-1.6935657) q[2];
sx q[2];
rz(-0.4099561) q[2];
rz(2.3112467) q[3];
sx q[3];
rz(-0.94873077) q[3];
sx q[3];
rz(1.4478987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53772563) q[0];
sx q[0];
rz(-0.69320885) q[0];
sx q[0];
rz(-2.895288) q[0];
rz(-2.314997) q[1];
sx q[1];
rz(-1.7431424) q[1];
sx q[1];
rz(2.8776317) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9519099) q[0];
sx q[0];
rz(-0.095746843) q[0];
sx q[0];
rz(1.4950947) q[0];
x q[1];
rz(2.6648681) q[2];
sx q[2];
rz(-0.84273224) q[2];
sx q[2];
rz(1.5155033) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8249656) q[1];
sx q[1];
rz(-2.3967603) q[1];
sx q[1];
rz(-1.1492351) q[1];
x q[2];
rz(-0.58707763) q[3];
sx q[3];
rz(-1.9403321) q[3];
sx q[3];
rz(2.0943506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0366514) q[2];
sx q[2];
rz(-1.2678601) q[2];
sx q[2];
rz(-0.71511739) q[2];
rz(-0.07130833) q[3];
sx q[3];
rz(-1.5569867) q[3];
sx q[3];
rz(2.3947072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0095373) q[0];
sx q[0];
rz(-1.4893463) q[0];
sx q[0];
rz(-2.706053) q[0];
rz(0.14777331) q[1];
sx q[1];
rz(-1.4316033) q[1];
sx q[1];
rz(-1.4685644) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4184121) q[0];
sx q[0];
rz(-1.3910798) q[0];
sx q[0];
rz(2.7030858) q[0];
rz(-1.4917489) q[2];
sx q[2];
rz(-2.2666396) q[2];
sx q[2];
rz(0.33338132) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.802001) q[1];
sx q[1];
rz(-2.7164384) q[1];
sx q[1];
rz(1.156686) q[1];
rz(2.7379509) q[3];
sx q[3];
rz(-1.3714681) q[3];
sx q[3];
rz(-2.2011873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.396951) q[2];
sx q[2];
rz(-1.7640742) q[2];
sx q[2];
rz(-0.89293876) q[2];
rz(1.2785771) q[3];
sx q[3];
rz(-1.9847001) q[3];
sx q[3];
rz(-1.6781835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29499149) q[0];
sx q[0];
rz(-2.8104267) q[0];
sx q[0];
rz(0.60605979) q[0];
rz(3.1312969) q[1];
sx q[1];
rz(-0.98774397) q[1];
sx q[1];
rz(0.079924718) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91658501) q[0];
sx q[0];
rz(-1.8226624) q[0];
sx q[0];
rz(-0.99117898) q[0];
rz(-pi) q[1];
rz(2.299231) q[2];
sx q[2];
rz(-1.3536705) q[2];
sx q[2];
rz(-2.6922243) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8680758) q[1];
sx q[1];
rz(-0.40765992) q[1];
sx q[1];
rz(0.056053921) q[1];
x q[2];
rz(1.3402088) q[3];
sx q[3];
rz(-1.7713431) q[3];
sx q[3];
rz(1.6236562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.14110485) q[2];
sx q[2];
rz(-2.2147369) q[2];
sx q[2];
rz(1.0151939) q[2];
rz(0.60123932) q[3];
sx q[3];
rz(-0.70178086) q[3];
sx q[3];
rz(-2.0950441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4973608) q[0];
sx q[0];
rz(-1.5174706) q[0];
sx q[0];
rz(0.68751412) q[0];
rz(0.58912206) q[1];
sx q[1];
rz(-0.67710572) q[1];
sx q[1];
rz(2.6893375) q[1];
rz(-1.4745787) q[2];
sx q[2];
rz(-2.0141891) q[2];
sx q[2];
rz(-0.63465848) q[2];
rz(-0.11832931) q[3];
sx q[3];
rz(-2.31349) q[3];
sx q[3];
rz(1.3009664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
