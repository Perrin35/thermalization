OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0107083) q[0];
sx q[0];
rz(3.8151725) q[0];
sx q[0];
rz(7.1060915) q[0];
rz(1.4505439) q[1];
sx q[1];
rz(-0.6757285) q[1];
sx q[1];
rz(-2.9339209) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0820706) q[0];
sx q[0];
rz(-0.86468177) q[0];
sx q[0];
rz(-0.45990277) q[0];
x q[1];
rz(2.7960289) q[2];
sx q[2];
rz(-1.0394382) q[2];
sx q[2];
rz(3.1329581) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.339141) q[1];
sx q[1];
rz(-1.4467518) q[1];
sx q[1];
rz(-1.5088827) q[1];
rz(-pi) q[2];
x q[2];
rz(0.071524044) q[3];
sx q[3];
rz(-1.2010964) q[3];
sx q[3];
rz(-0.92310753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9660008) q[2];
sx q[2];
rz(-1.5755743) q[2];
sx q[2];
rz(1.2292181) q[2];
rz(-1.843533) q[3];
sx q[3];
rz(-1.3763873) q[3];
sx q[3];
rz(-1.957533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9640279) q[0];
sx q[0];
rz(-0.25968817) q[0];
sx q[0];
rz(-0.92726707) q[0];
rz(-0.19993965) q[1];
sx q[1];
rz(-1.4727458) q[1];
sx q[1];
rz(2.1520069) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9703116) q[0];
sx q[0];
rz(-1.6785445) q[0];
sx q[0];
rz(1.9338444) q[0];
x q[1];
rz(1.5275498) q[2];
sx q[2];
rz(-2.0048365) q[2];
sx q[2];
rz(-2.1021646) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.67827077) q[1];
sx q[1];
rz(-0.090677977) q[1];
sx q[1];
rz(-2.1342127) q[1];
rz(-2.5603676) q[3];
sx q[3];
rz(-2.1970501) q[3];
sx q[3];
rz(-3.0581829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2135311) q[2];
sx q[2];
rz(-1.7240883) q[2];
sx q[2];
rz(0.35378635) q[2];
rz(2.1900322) q[3];
sx q[3];
rz(-2.3944201) q[3];
sx q[3];
rz(0.32683867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3554409) q[0];
sx q[0];
rz(-0.94796258) q[0];
sx q[0];
rz(2.1308664) q[0];
rz(-3.082869) q[1];
sx q[1];
rz(-1.6207691) q[1];
sx q[1];
rz(-2.7322863) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.302251) q[0];
sx q[0];
rz(-0.88345655) q[0];
sx q[0];
rz(0.99498596) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9818166) q[2];
sx q[2];
rz(-1.7226698) q[2];
sx q[2];
rz(1.077026) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9524556) q[1];
sx q[1];
rz(-0.5824832) q[1];
sx q[1];
rz(0.36242318) q[1];
x q[2];
rz(0.50893754) q[3];
sx q[3];
rz(-1.5526315) q[3];
sx q[3];
rz(-3.09336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0851486) q[2];
sx q[2];
rz(-1.2480382) q[2];
sx q[2];
rz(0.15963456) q[2];
rz(1.8917482) q[3];
sx q[3];
rz(-1.4188473) q[3];
sx q[3];
rz(-2.0554525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98642629) q[0];
sx q[0];
rz(-0.43743375) q[0];
sx q[0];
rz(1.0945818) q[0];
rz(-1.9704341) q[1];
sx q[1];
rz(-1.1181701) q[1];
sx q[1];
rz(2.1280682) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5399234) q[0];
sx q[0];
rz(-0.058383103) q[0];
sx q[0];
rz(-2.4235382) q[0];
rz(-0.0036448467) q[2];
sx q[2];
rz(-2.0652152) q[2];
sx q[2];
rz(1.0629176) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3588874) q[1];
sx q[1];
rz(-0.26391477) q[1];
sx q[1];
rz(0.15586075) q[1];
rz(-1.1356527) q[3];
sx q[3];
rz(-1.6683104) q[3];
sx q[3];
rz(-1.9103736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1363498) q[2];
sx q[2];
rz(-1.3815657) q[2];
sx q[2];
rz(-1.8278149) q[2];
rz(2.2037196) q[3];
sx q[3];
rz(-1.8504986) q[3];
sx q[3];
rz(-0.098793678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90792847) q[0];
sx q[0];
rz(-1.6133244) q[0];
sx q[0];
rz(2.0966356) q[0];
rz(-0.16432556) q[1];
sx q[1];
rz(-1.798809) q[1];
sx q[1];
rz(-0.16990653) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63498492) q[0];
sx q[0];
rz(-1.7568551) q[0];
sx q[0];
rz(3.0691514) q[0];
x q[1];
rz(0.81973524) q[2];
sx q[2];
rz(-2.6728874) q[2];
sx q[2];
rz(-0.021989487) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.85346881) q[1];
sx q[1];
rz(-2.2125724) q[1];
sx q[1];
rz(0.30499129) q[1];
rz(0.72615926) q[3];
sx q[3];
rz(-1.4724331) q[3];
sx q[3];
rz(-1.5226328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5728411) q[2];
sx q[2];
rz(-0.70256394) q[2];
sx q[2];
rz(2.11002) q[2];
rz(-2.0555563) q[3];
sx q[3];
rz(-2.3417754) q[3];
sx q[3];
rz(-1.4097376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3373435) q[0];
sx q[0];
rz(-2.0954837) q[0];
sx q[0];
rz(1.7403437) q[0];
rz(-2.7119472) q[1];
sx q[1];
rz(-0.88637561) q[1];
sx q[1];
rz(-1.8362129) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1032377) q[0];
sx q[0];
rz(-1.2631386) q[0];
sx q[0];
rz(-2.3820113) q[0];
rz(3.1370509) q[2];
sx q[2];
rz(-2.0097187) q[2];
sx q[2];
rz(2.2786841) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.77338058) q[1];
sx q[1];
rz(-1.6016593) q[1];
sx q[1];
rz(-2.8931151) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8135728) q[3];
sx q[3];
rz(-1.7194028) q[3];
sx q[3];
rz(0.44156238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1244916) q[2];
sx q[2];
rz(-1.2064826) q[2];
sx q[2];
rz(-2.1567832) q[2];
rz(1.5752327) q[3];
sx q[3];
rz(-1.5519578) q[3];
sx q[3];
rz(2.8563833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8816836) q[0];
sx q[0];
rz(-2.8398828) q[0];
sx q[0];
rz(1.786422) q[0];
rz(-3.1175218) q[1];
sx q[1];
rz(-1.5211952) q[1];
sx q[1];
rz(-0.37758652) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1771072) q[0];
sx q[0];
rz(-1.4146574) q[0];
sx q[0];
rz(0.6556169) q[0];
rz(2.2607498) q[2];
sx q[2];
rz(-2.8683148) q[2];
sx q[2];
rz(-0.6169332) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7086975) q[1];
sx q[1];
rz(-2.5835977) q[1];
sx q[1];
rz(-0.79496986) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0041647) q[3];
sx q[3];
rz(-1.2655756) q[3];
sx q[3];
rz(0.63852541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5923803) q[2];
sx q[2];
rz(-0.33752957) q[2];
sx q[2];
rz(0.40360061) q[2];
rz(-2.9523383) q[3];
sx q[3];
rz(-2.0892102) q[3];
sx q[3];
rz(2.6846867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52983487) q[0];
sx q[0];
rz(-2.5921322) q[0];
sx q[0];
rz(-2.8305565) q[0];
rz(0.95651904) q[1];
sx q[1];
rz(-1.4776968) q[1];
sx q[1];
rz(2.8172353) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7081147) q[0];
sx q[0];
rz(-1.38033) q[0];
sx q[0];
rz(-2.6300927) q[0];
rz(-0.51028549) q[2];
sx q[2];
rz(-1.6225177) q[2];
sx q[2];
rz(3.1003458) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.16087577) q[1];
sx q[1];
rz(-2.389713) q[1];
sx q[1];
rz(-2.1999825) q[1];
rz(-pi) q[2];
rz(-1.1127541) q[3];
sx q[3];
rz(-2.4107217) q[3];
sx q[3];
rz(2.3422086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4057464) q[2];
sx q[2];
rz(-2.3509071) q[2];
sx q[2];
rz(2.9131367) q[2];
rz(-0.53260803) q[3];
sx q[3];
rz(-2.3753128) q[3];
sx q[3];
rz(0.84958357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.547895) q[0];
sx q[0];
rz(-0.51012796) q[0];
sx q[0];
rz(1.8713895) q[0];
rz(-0.58865976) q[1];
sx q[1];
rz(-1.9344067) q[1];
sx q[1];
rz(0.38280815) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24341881) q[0];
sx q[0];
rz(-2.9110258) q[0];
sx q[0];
rz(-1.4617993) q[0];
rz(-pi) q[1];
rz(2.7076376) q[2];
sx q[2];
rz(-1.343691) q[2];
sx q[2];
rz(1.8411098) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8296903) q[1];
sx q[1];
rz(-1.4451348) q[1];
sx q[1];
rz(-2.4290393) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0625913) q[3];
sx q[3];
rz(-1.8100836) q[3];
sx q[3];
rz(-0.21904473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1560912) q[2];
sx q[2];
rz(-1.5727377) q[2];
sx q[2];
rz(2.7413979) q[2];
rz(0.24724809) q[3];
sx q[3];
rz(-1.4618382) q[3];
sx q[3];
rz(-0.39321536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3006712) q[0];
sx q[0];
rz(-1.8074169) q[0];
sx q[0];
rz(-1.0155431) q[0];
rz(-0.66889846) q[1];
sx q[1];
rz(-1.4543507) q[1];
sx q[1];
rz(-1.6355754) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72476232) q[0];
sx q[0];
rz(-2.2106631) q[0];
sx q[0];
rz(-0.93389966) q[0];
rz(2.4898275) q[2];
sx q[2];
rz(-0.6547857) q[2];
sx q[2];
rz(0.98818615) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2336894) q[1];
sx q[1];
rz(-2.0504867) q[1];
sx q[1];
rz(-3.0976899) q[1];
rz(-pi) q[2];
rz(-0.90088441) q[3];
sx q[3];
rz(-0.79798079) q[3];
sx q[3];
rz(2.1547998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9754535) q[2];
sx q[2];
rz(-1.0756805) q[2];
sx q[2];
rz(-1.5258741) q[2];
rz(0.29159355) q[3];
sx q[3];
rz(-0.44277954) q[3];
sx q[3];
rz(0.35167882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46353729) q[0];
sx q[0];
rz(-0.96374496) q[0];
sx q[0];
rz(-2.5441334) q[0];
rz(2.8529104) q[1];
sx q[1];
rz(-2.3434227) q[1];
sx q[1];
rz(-3.1176288) q[1];
rz(2.504442) q[2];
sx q[2];
rz(-0.49497866) q[2];
sx q[2];
rz(-0.53866932) q[2];
rz(0.45997672) q[3];
sx q[3];
rz(-1.796985) q[3];
sx q[3];
rz(-1.7225869) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
