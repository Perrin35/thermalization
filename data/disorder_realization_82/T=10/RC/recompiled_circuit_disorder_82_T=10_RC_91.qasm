OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80469552) q[0];
sx q[0];
rz(-1.0372294) q[0];
sx q[0];
rz(-2.7859935) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(-1.27966) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1188388) q[0];
sx q[0];
rz(-1.4316598) q[0];
sx q[0];
rz(-1.3214146) q[0];
rz(-pi) q[1];
rz(-1.7945292) q[2];
sx q[2];
rz(-1.1587843) q[2];
sx q[2];
rz(1.7681233) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3763334) q[1];
sx q[1];
rz(-1.9783124) q[1];
sx q[1];
rz(2.5024662) q[1];
rz(2.9079307) q[3];
sx q[3];
rz(-1.3705001) q[3];
sx q[3];
rz(-1.3113126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1352284) q[2];
sx q[2];
rz(-0.93868119) q[2];
sx q[2];
rz(0.65650666) q[2];
rz(-2.4025829) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(-2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1332557) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(2.546229) q[0];
rz(-3.0796675) q[1];
sx q[1];
rz(-1.9134816) q[1];
sx q[1];
rz(0.48746902) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2330403) q[0];
sx q[0];
rz(-1.6516343) q[0];
sx q[0];
rz(3.0768865) q[0];
x q[1];
rz(2.3357453) q[2];
sx q[2];
rz(-0.57493756) q[2];
sx q[2];
rz(0.88428674) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4844955) q[1];
sx q[1];
rz(-2.8946981) q[1];
sx q[1];
rz(-0.8685649) q[1];
rz(-pi) q[2];
rz(-1.3120193) q[3];
sx q[3];
rz(-2.8142455) q[3];
sx q[3];
rz(2.5395288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9618824) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(-0.3953735) q[2];
rz(-2.0987434) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(0.038671967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19668002) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(-0.19038598) q[0];
rz(-0.12292513) q[1];
sx q[1];
rz(-0.38750896) q[1];
sx q[1];
rz(2.9188459) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73244625) q[0];
sx q[0];
rz(-1.961686) q[0];
sx q[0];
rz(-1.976165) q[0];
x q[1];
rz(2.9163755) q[2];
sx q[2];
rz(-1.2080964) q[2];
sx q[2];
rz(-2.1622216) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0157156) q[1];
sx q[1];
rz(-1.5029969) q[1];
sx q[1];
rz(2.5281639) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1915216) q[3];
sx q[3];
rz(-0.40682236) q[3];
sx q[3];
rz(-2.1086958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8824076) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(-1.9474691) q[2];
rz(-1.0549226) q[3];
sx q[3];
rz(-1.4849562) q[3];
sx q[3];
rz(-2.1779493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7524183) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(-1.4032723) q[0];
rz(1.6482884) q[1];
sx q[1];
rz(-1.5999258) q[1];
sx q[1];
rz(-0.9202252) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.614914) q[0];
sx q[0];
rz(-0.81256142) q[0];
sx q[0];
rz(2.148669) q[0];
rz(-pi) q[1];
rz(-2.9158981) q[2];
sx q[2];
rz(-2.7719471) q[2];
sx q[2];
rz(-2.8452578) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.11776609) q[1];
sx q[1];
rz(-0.95925602) q[1];
sx q[1];
rz(2.141342) q[1];
rz(-1.7761049) q[3];
sx q[3];
rz(-1.218154) q[3];
sx q[3];
rz(-2.7384788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.197864) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(-1.2801923) q[2];
rz(-2.3222893) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.458805) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(-1.4047594) q[0];
rz(-2.3732896) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(-0.4531025) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8720854) q[0];
sx q[0];
rz(-0.082490248) q[0];
sx q[0];
rz(2.0399658) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4278973) q[2];
sx q[2];
rz(-1.5362527) q[2];
sx q[2];
rz(-2.3847716) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.52757712) q[1];
sx q[1];
rz(-1.3506538) q[1];
sx q[1];
rz(3.0575391) q[1];
rz(0.94282486) q[3];
sx q[3];
rz(-0.70126611) q[3];
sx q[3];
rz(-3.0758465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.12864628) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(-1.6298693) q[2];
rz(2.417918) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(-2.2912912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1081651) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(0.23813716) q[0];
rz(-2.761633) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(1.5135117) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0673163) q[0];
sx q[0];
rz(-2.5551676) q[0];
sx q[0];
rz(2.5308454) q[0];
rz(-pi) q[1];
rz(2.020535) q[2];
sx q[2];
rz(-1.5151086) q[2];
sx q[2];
rz(-2.4186717) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7102393) q[1];
sx q[1];
rz(-0.87190404) q[1];
sx q[1];
rz(0.61183521) q[1];
x q[2];
rz(2.2699039) q[3];
sx q[3];
rz(-2.5786434) q[3];
sx q[3];
rz(-0.20760078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0825519) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(1.9980105) q[2];
rz(-0.13051662) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1259574) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(-2.7440199) q[0];
rz(-1.3145087) q[1];
sx q[1];
rz(-2.59771) q[1];
sx q[1];
rz(-1.9546753) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8531187) q[0];
sx q[0];
rz(-1.0736199) q[0];
sx q[0];
rz(-2.8209646) q[0];
x q[1];
rz(-1.5415807) q[2];
sx q[2];
rz(-1.6299106) q[2];
sx q[2];
rz(0.37994775) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.742332) q[1];
sx q[1];
rz(-2.2982344) q[1];
sx q[1];
rz(-2.9620693) q[1];
rz(-pi) q[2];
rz(2.5116634) q[3];
sx q[3];
rz(-1.3327206) q[3];
sx q[3];
rz(1.5743953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26178965) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(-1.3558033) q[2];
rz(1.6342182) q[3];
sx q[3];
rz(-2.5528788) q[3];
sx q[3];
rz(-0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(2.3223406) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(-0.85574714) q[0];
rz(0.021727173) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(-1.0303248) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4521342) q[0];
sx q[0];
rz(-1.2297213) q[0];
sx q[0];
rz(2.336691) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4872929) q[2];
sx q[2];
rz(-0.52831542) q[2];
sx q[2];
rz(-2.5022142) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9955666) q[1];
sx q[1];
rz(-1.4194173) q[1];
sx q[1];
rz(-1.2703018) q[1];
rz(-pi) q[2];
rz(0.86352591) q[3];
sx q[3];
rz(-2.415495) q[3];
sx q[3];
rz(-1.8973779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8490303) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(3.0714152) q[2];
rz(2.3146546) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(-0.58644811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4940015) q[0];
sx q[0];
rz(-1.4275455) q[0];
sx q[0];
rz(2.8253187) q[0];
rz(2.0896185) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(-2.0064328) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59199698) q[0];
sx q[0];
rz(-1.9599008) q[0];
sx q[0];
rz(-1.367021) q[0];
rz(-pi) q[1];
rz(1.4831545) q[2];
sx q[2];
rz(-0.32155514) q[2];
sx q[2];
rz(1.9784387) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.89716298) q[1];
sx q[1];
rz(-2.0532554) q[1];
sx q[1];
rz(-0.06628118) q[1];
rz(-pi) q[2];
rz(-3.1177403) q[3];
sx q[3];
rz(-2.1937222) q[3];
sx q[3];
rz(0.21104392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6241374) q[2];
sx q[2];
rz(-2.3949261) q[2];
sx q[2];
rz(-1.5967782) q[2];
rz(-2.4638713) q[3];
sx q[3];
rz(-0.8422519) q[3];
sx q[3];
rz(-2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90821663) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(-2.7897575) q[0];
rz(-2.8219163) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(0.19616729) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10712121) q[0];
sx q[0];
rz(-2.2959318) q[0];
sx q[0];
rz(2.3496773) q[0];
rz(-1.9253795) q[2];
sx q[2];
rz(-2.9973929) q[2];
sx q[2];
rz(1.2186288) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0077121) q[1];
sx q[1];
rz(-0.73426437) q[1];
sx q[1];
rz(-1.5478565) q[1];
x q[2];
rz(-0.84368002) q[3];
sx q[3];
rz(-2.0920678) q[3];
sx q[3];
rz(-2.664444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3907884) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(0.081136726) q[2];
rz(1.9598512) q[3];
sx q[3];
rz(-2.2406082) q[3];
sx q[3];
rz(-1.0749764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.022973013) q[0];
sx q[0];
rz(-1.5587627) q[0];
sx q[0];
rz(-1.8297304) q[0];
rz(1.8267869) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(1.6191471) q[2];
sx q[2];
rz(-1.3370677) q[2];
sx q[2];
rz(1.7481902) q[2];
rz(-1.374791) q[3];
sx q[3];
rz(-0.23402611) q[3];
sx q[3];
rz(0.2454161) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
