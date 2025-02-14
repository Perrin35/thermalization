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
rz(1.243408) q[0];
sx q[0];
rz(4.7060407) q[0];
sx q[0];
rz(10.336182) q[0];
rz(2.4721594) q[1];
sx q[1];
rz(3.4179847) q[1];
sx q[1];
rz(8.5951947) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88032915) q[0];
sx q[0];
rz(-1.1733049) q[0];
sx q[0];
rz(0.3573768) q[0];
rz(-pi) q[1];
rz(-2.787399) q[2];
sx q[2];
rz(-2.8787432) q[2];
sx q[2];
rz(2.7517002) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3717229) q[1];
sx q[1];
rz(-0.81156213) q[1];
sx q[1];
rz(-2.7219591) q[1];
x q[2];
rz(-0.43609332) q[3];
sx q[3];
rz(-1.4319518) q[3];
sx q[3];
rz(0.88263765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7726045) q[2];
sx q[2];
rz(-2.172894) q[2];
sx q[2];
rz(-1.1645092) q[2];
rz(0.72689593) q[3];
sx q[3];
rz(-1.8098857) q[3];
sx q[3];
rz(-0.070778457) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1233391) q[0];
sx q[0];
rz(-2.4788719) q[0];
sx q[0];
rz(2.8501999) q[0];
rz(0.60028752) q[1];
sx q[1];
rz(-0.95708668) q[1];
sx q[1];
rz(-0.16500638) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5614189) q[0];
sx q[0];
rz(-0.2720755) q[0];
sx q[0];
rz(-1.838057) q[0];
rz(-pi) q[1];
rz(1.3403471) q[2];
sx q[2];
rz(-2.682095) q[2];
sx q[2];
rz(0.27333096) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.2272039) q[1];
sx q[1];
rz(-1.640135) q[1];
sx q[1];
rz(-0.68080618) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6253774) q[3];
sx q[3];
rz(-1.8187869) q[3];
sx q[3];
rz(2.8021106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0642285) q[2];
sx q[2];
rz(-2.1891429) q[2];
sx q[2];
rz(-2.9020818) q[2];
rz(0.32660487) q[3];
sx q[3];
rz(-2.9686847) q[3];
sx q[3];
rz(-3.0467564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8851682) q[0];
sx q[0];
rz(-0.435985) q[0];
sx q[0];
rz(-1.5727795) q[0];
rz(0.39372152) q[1];
sx q[1];
rz(-0.55501333) q[1];
sx q[1];
rz(-0.47600019) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25903738) q[0];
sx q[0];
rz(-0.60094423) q[0];
sx q[0];
rz(-1.9447692) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4127998) q[2];
sx q[2];
rz(-2.9523628) q[2];
sx q[2];
rz(0.10052027) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.90895528) q[1];
sx q[1];
rz(-0.74027762) q[1];
sx q[1];
rz(2.2038682) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2331486) q[3];
sx q[3];
rz(-1.3364571) q[3];
sx q[3];
rz(0.657814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8126882) q[2];
sx q[2];
rz(-0.71949553) q[2];
sx q[2];
rz(-2.1997814) q[2];
rz(-1.4224667) q[3];
sx q[3];
rz(-2.7673281) q[3];
sx q[3];
rz(-0.67371887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0734237) q[0];
sx q[0];
rz(-2.9015559) q[0];
sx q[0];
rz(1.5330676) q[0];
rz(0.34670058) q[1];
sx q[1];
rz(-2.0997212) q[1];
sx q[1];
rz(2.9516721) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8297716) q[0];
sx q[0];
rz(-1.8320727) q[0];
sx q[0];
rz(2.3992062) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0713351) q[2];
sx q[2];
rz(-0.53218666) q[2];
sx q[2];
rz(-2.4376873) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0611113) q[1];
sx q[1];
rz(-0.60855344) q[1];
sx q[1];
rz(-1.4707277) q[1];
rz(-pi) q[2];
rz(-1.1433825) q[3];
sx q[3];
rz(-1.1955002) q[3];
sx q[3];
rz(2.0209598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7803663) q[2];
sx q[2];
rz(-2.2915338) q[2];
sx q[2];
rz(0.38193199) q[2];
rz(0.10968883) q[3];
sx q[3];
rz(-1.0332801) q[3];
sx q[3];
rz(-0.1194574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0383976) q[0];
sx q[0];
rz(-1.1829475) q[0];
sx q[0];
rz(2.2947445) q[0];
rz(3.0021744) q[1];
sx q[1];
rz(-0.81652313) q[1];
sx q[1];
rz(-0.74904186) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3232305) q[0];
sx q[0];
rz(-1.4233755) q[0];
sx q[0];
rz(2.4458103) q[0];
rz(-pi) q[1];
rz(0.26450677) q[2];
sx q[2];
rz(-0.49969426) q[2];
sx q[2];
rz(-3.0124912) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.79758039) q[1];
sx q[1];
rz(-1.3770119) q[1];
sx q[1];
rz(1.4414201) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0723125) q[3];
sx q[3];
rz(-1.673888) q[3];
sx q[3];
rz(0.13971288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.9355115) q[2];
sx q[2];
rz(-1.2568018) q[2];
sx q[2];
rz(2.0743267) q[2];
rz(0.71077985) q[3];
sx q[3];
rz(-1.4830517) q[3];
sx q[3];
rz(1.7162286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34273219) q[0];
sx q[0];
rz(-1.6628168) q[0];
sx q[0];
rz(2.1730098) q[0];
rz(-0.69106483) q[1];
sx q[1];
rz(-1.7993118) q[1];
sx q[1];
rz(-1.3074646) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1145612) q[0];
sx q[0];
rz(-1.3163211) q[0];
sx q[0];
rz(2.4976394) q[0];
rz(-pi) q[1];
rz(0.5796104) q[2];
sx q[2];
rz(-0.34380772) q[2];
sx q[2];
rz(-1.0506309) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4279514) q[1];
sx q[1];
rz(-0.390807) q[1];
sx q[1];
rz(1.7094088) q[1];
x q[2];
rz(2.9581466) q[3];
sx q[3];
rz(-1.8357092) q[3];
sx q[3];
rz(-1.9062756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.4751733) q[2];
sx q[2];
rz(-0.23491493) q[2];
sx q[2];
rz(1.2972181) q[2];
rz(1.3951067) q[3];
sx q[3];
rz(-1.3336983) q[3];
sx q[3];
rz(-1.6102128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7995826) q[0];
sx q[0];
rz(-0.81556773) q[0];
sx q[0];
rz(-0.94014257) q[0];
rz(-2.4932585) q[1];
sx q[1];
rz(-1.9377361) q[1];
sx q[1];
rz(-0.26888332) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7358582) q[0];
sx q[0];
rz(-2.5395505) q[0];
sx q[0];
rz(-2.6989537) q[0];
rz(-pi) q[1];
rz(-2.778947) q[2];
sx q[2];
rz(-0.73604903) q[2];
sx q[2];
rz(-0.48378746) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2317048) q[1];
sx q[1];
rz(-2.1429166) q[1];
sx q[1];
rz(2.8673299) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.712065) q[3];
sx q[3];
rz(-1.2614354) q[3];
sx q[3];
rz(-1.3573028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.0073283422) q[2];
sx q[2];
rz(-1.8625883) q[2];
sx q[2];
rz(-2.4556665) q[2];
rz(2.4053597) q[3];
sx q[3];
rz(-2.1015621) q[3];
sx q[3];
rz(0.22549103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69970423) q[0];
sx q[0];
rz(-1.2208953) q[0];
sx q[0];
rz(0.7487444) q[0];
rz(3.0486095) q[1];
sx q[1];
rz(-2.3817606) q[1];
sx q[1];
rz(1.071113) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1919928) q[0];
sx q[0];
rz(-1.1566391) q[0];
sx q[0];
rz(0.55121514) q[0];
rz(-pi) q[1];
rz(0.65437859) q[2];
sx q[2];
rz(-1.9438231) q[2];
sx q[2];
rz(-1.6668863) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.86102137) q[1];
sx q[1];
rz(-1.264466) q[1];
sx q[1];
rz(0.88545274) q[1];
x q[2];
rz(1.6900916) q[3];
sx q[3];
rz(-0.41364851) q[3];
sx q[3];
rz(-0.0095490015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6894655) q[2];
sx q[2];
rz(-2.5048544) q[2];
sx q[2];
rz(-0.20991906) q[2];
rz(-0.12957761) q[3];
sx q[3];
rz(-2.1942873) q[3];
sx q[3];
rz(1.5773704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.121948) q[0];
sx q[0];
rz(-2.8804998) q[0];
sx q[0];
rz(-3.1280532) q[0];
rz(-0.53567046) q[1];
sx q[1];
rz(-0.56709138) q[1];
sx q[1];
rz(3.1279235) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5608198) q[0];
sx q[0];
rz(-0.81308621) q[0];
sx q[0];
rz(-1.2570501) q[0];
rz(-pi) q[1];
rz(-1.0837567) q[2];
sx q[2];
rz(-1.221861) q[2];
sx q[2];
rz(-1.2020122) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7175258) q[1];
sx q[1];
rz(-0.88554157) q[1];
sx q[1];
rz(1.9937033) q[1];
rz(-pi) q[2];
rz(2.3914371) q[3];
sx q[3];
rz(-1.8492846) q[3];
sx q[3];
rz(0.089547308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.51836625) q[2];
sx q[2];
rz(-1.7609111) q[2];
sx q[2];
rz(-1.3631442) q[2];
rz(-0.81950435) q[3];
sx q[3];
rz(-2.6624811) q[3];
sx q[3];
rz(-0.68377408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0226456) q[0];
sx q[0];
rz(-2.2830257) q[0];
sx q[0];
rz(1.488142) q[0];
rz(0.34291521) q[1];
sx q[1];
rz(-1.5838793) q[1];
sx q[1];
rz(0.75103474) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9404749) q[0];
sx q[0];
rz(-1.6342666) q[0];
sx q[0];
rz(1.2918143) q[0];
rz(2.033634) q[2];
sx q[2];
rz(-1.3909391) q[2];
sx q[2];
rz(-0.14000574) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18606762) q[1];
sx q[1];
rz(-2.1677818) q[1];
sx q[1];
rz(1.7138193) q[1];
rz(-pi) q[2];
rz(-1.6609825) q[3];
sx q[3];
rz(-1.6038365) q[3];
sx q[3];
rz(0.15109135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47120961) q[2];
sx q[2];
rz(-0.94380108) q[2];
sx q[2];
rz(0.38718265) q[2];
rz(0.47532982) q[3];
sx q[3];
rz(-1.9742222) q[3];
sx q[3];
rz(-2.3502684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8925856) q[0];
sx q[0];
rz(-0.96207608) q[0];
sx q[0];
rz(0.67208653) q[0];
rz(-2.7723906) q[1];
sx q[1];
rz(-0.75405706) q[1];
sx q[1];
rz(-1.3934607) q[1];
rz(-1.3849707) q[2];
sx q[2];
rz(-0.86238774) q[2];
sx q[2];
rz(2.9288616) q[2];
rz(-3.0283916) q[3];
sx q[3];
rz(-1.8424215) q[3];
sx q[3];
rz(2.7383534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
