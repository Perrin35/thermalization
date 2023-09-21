OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6931273) q[0];
sx q[0];
rz(-0.52283302) q[0];
sx q[0];
rz(-0.62358207) q[0];
rz(-2.8514255) q[1];
sx q[1];
rz(-0.71915141) q[1];
sx q[1];
rz(-2.6410988) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96007632) q[0];
sx q[0];
rz(-1.5746563) q[0];
sx q[0];
rz(-0.5998248) q[0];
x q[1];
rz(1.7206465) q[2];
sx q[2];
rz(-1.4784001) q[2];
sx q[2];
rz(2.7211651) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7841255) q[1];
sx q[1];
rz(-0.81175121) q[1];
sx q[1];
rz(-1.1151421) q[1];
rz(-0.15668232) q[3];
sx q[3];
rz(-2.0041668) q[3];
sx q[3];
rz(2.3436848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3502675) q[2];
sx q[2];
rz(-1.9359549) q[2];
sx q[2];
rz(-1.9187437) q[2];
rz(1.4482927) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(-2.1712415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7704849) q[0];
sx q[0];
rz(-1.4828232) q[0];
sx q[0];
rz(-2.0626542) q[0];
rz(1.3868015) q[1];
sx q[1];
rz(-0.81258041) q[1];
sx q[1];
rz(0.66545495) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5011713) q[0];
sx q[0];
rz(-1.8542395) q[0];
sx q[0];
rz(-2.6675176) q[0];
rz(-pi) q[1];
rz(2.6639054) q[2];
sx q[2];
rz(-2.7825232) q[2];
sx q[2];
rz(1.2815086) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.39311073) q[1];
sx q[1];
rz(-2.0512274) q[1];
sx q[1];
rz(1.837681) q[1];
rz(-pi) q[2];
rz(-0.86136787) q[3];
sx q[3];
rz(-1.5966166) q[3];
sx q[3];
rz(-2.4840419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5370496) q[2];
sx q[2];
rz(-1.4031354) q[2];
sx q[2];
rz(0.075142168) q[2];
rz(1.4705307) q[3];
sx q[3];
rz(-1.1276779) q[3];
sx q[3];
rz(-2.4501734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5866518) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(0.75876045) q[0];
rz(1.2930019) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(-1.4000777) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6436359) q[0];
sx q[0];
rz(-2.6640577) q[0];
sx q[0];
rz(0.60995539) q[0];
rz(-pi) q[1];
rz(1.4351575) q[2];
sx q[2];
rz(-1.3856158) q[2];
sx q[2];
rz(0.26091012) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2565508) q[1];
sx q[1];
rz(-2.1532144) q[1];
sx q[1];
rz(0.38862733) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2657884) q[3];
sx q[3];
rz(-1.4121571) q[3];
sx q[3];
rz(0.97914417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.65163461) q[2];
sx q[2];
rz(-2.659446) q[2];
sx q[2];
rz(0.65762323) q[2];
rz(1.1714606) q[3];
sx q[3];
rz(-1.6572584) q[3];
sx q[3];
rz(1.4191779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79384971) q[0];
sx q[0];
rz(-2.1934953) q[0];
sx q[0];
rz(-1.6963652) q[0];
rz(-1.6943278) q[1];
sx q[1];
rz(-1.4935962) q[1];
sx q[1];
rz(-0.34805527) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6570243) q[0];
sx q[0];
rz(-2.1977402) q[0];
sx q[0];
rz(-1.2059962) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15050998) q[2];
sx q[2];
rz(-1.2005271) q[2];
sx q[2];
rz(2.9589257) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0749803) q[1];
sx q[1];
rz(-2.1790824) q[1];
sx q[1];
rz(0.10886701) q[1];
rz(-pi) q[2];
rz(2.4481941) q[3];
sx q[3];
rz(-2.9326673) q[3];
sx q[3];
rz(0.25017504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7248914) q[2];
sx q[2];
rz(-1.3092224) q[2];
sx q[2];
rz(-2.7187738) q[2];
rz(0.73741284) q[3];
sx q[3];
rz(-2.335572) q[3];
sx q[3];
rz(3.056934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87930644) q[0];
sx q[0];
rz(-1.3129741) q[0];
sx q[0];
rz(-0.53043956) q[0];
rz(-2.2166705) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(1.8431429) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86622483) q[0];
sx q[0];
rz(-0.21731649) q[0];
sx q[0];
rz(-2.2989681) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41579397) q[2];
sx q[2];
rz(-0.64877629) q[2];
sx q[2];
rz(-0.4817889) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7057719) q[1];
sx q[1];
rz(-2.768369) q[1];
sx q[1];
rz(2.865764) q[1];
x q[2];
rz(-0.97814822) q[3];
sx q[3];
rz(-1.1453298) q[3];
sx q[3];
rz(3.0654207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2003145) q[2];
sx q[2];
rz(-1.986074) q[2];
sx q[2];
rz(-0.47362622) q[2];
rz(0.099362699) q[3];
sx q[3];
rz(-1.8615581) q[3];
sx q[3];
rz(-2.30106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50399238) q[0];
sx q[0];
rz(-1.7853328) q[0];
sx q[0];
rz(2.1437058) q[0];
rz(2.2672794) q[1];
sx q[1];
rz(-1.0214146) q[1];
sx q[1];
rz(0.46674892) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40690639) q[0];
sx q[0];
rz(-0.91306251) q[0];
sx q[0];
rz(1.8440194) q[0];
rz(-pi) q[1];
rz(2.6830707) q[2];
sx q[2];
rz(-2.6715171) q[2];
sx q[2];
rz(0.93570659) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1280061) q[1];
sx q[1];
rz(-0.15833536) q[1];
sx q[1];
rz(-0.30551417) q[1];
rz(1.7344597) q[3];
sx q[3];
rz(-1.6657077) q[3];
sx q[3];
rz(-1.519219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2816887) q[2];
sx q[2];
rz(-0.47912654) q[2];
sx q[2];
rz(1.5647282) q[2];
rz(-0.62670296) q[3];
sx q[3];
rz(-1.7765216) q[3];
sx q[3];
rz(-2.4600162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8108114) q[0];
sx q[0];
rz(-2.2450876) q[0];
sx q[0];
rz(-2.7217641) q[0];
rz(2.9201674) q[1];
sx q[1];
rz(-2.6629993) q[1];
sx q[1];
rz(-1.5931169) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89824642) q[0];
sx q[0];
rz(-1.4833437) q[0];
sx q[0];
rz(1.5623708) q[0];
rz(2.1599342) q[2];
sx q[2];
rz(-0.33827153) q[2];
sx q[2];
rz(-0.32064082) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6429813) q[1];
sx q[1];
rz(-1.4652518) q[1];
sx q[1];
rz(-1.5555744) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7967954) q[3];
sx q[3];
rz(-2.0119917) q[3];
sx q[3];
rz(1.3519209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5856813) q[2];
sx q[2];
rz(-2.6101117) q[2];
sx q[2];
rz(-1.7377724) q[2];
rz(-2.3445271) q[3];
sx q[3];
rz(-2.6795487) q[3];
sx q[3];
rz(-1.2169303) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4306915) q[0];
sx q[0];
rz(-2.4276908) q[0];
sx q[0];
rz(-2.8523493) q[0];
rz(2.5166683) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(-1.3141059) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4411366) q[0];
sx q[0];
rz(-1.6807798) q[0];
sx q[0];
rz(2.3814047) q[0];
rz(-pi) q[1];
rz(-0.35347519) q[2];
sx q[2];
rz(-2.3620053) q[2];
sx q[2];
rz(-2.0445063) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5319547) q[1];
sx q[1];
rz(-1.6913809) q[1];
sx q[1];
rz(-3.0871255) q[1];
x q[2];
rz(0.81359158) q[3];
sx q[3];
rz(-1.3435257) q[3];
sx q[3];
rz(3.1384625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.73359314) q[2];
sx q[2];
rz(-1.5780129) q[2];
sx q[2];
rz(1.7129664) q[2];
rz(-0.96380487) q[3];
sx q[3];
rz(-1.0771841) q[3];
sx q[3];
rz(-1.0296286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52255094) q[0];
sx q[0];
rz(-1.4551117) q[0];
sx q[0];
rz(1.2458941) q[0];
rz(0.11101162) q[1];
sx q[1];
rz(-1.1975892) q[1];
sx q[1];
rz(0.54661173) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3205991) q[0];
sx q[0];
rz(-0.59251596) q[0];
sx q[0];
rz(2.2792363) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75066363) q[2];
sx q[2];
rz(-2.6235136) q[2];
sx q[2];
rz(0.93875611) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6555772) q[1];
sx q[1];
rz(-1.7421107) q[1];
sx q[1];
rz(-1.1737215) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6528477) q[3];
sx q[3];
rz(-2.1663323) q[3];
sx q[3];
rz(-2.3084156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1116011) q[2];
sx q[2];
rz(-2.2932055) q[2];
sx q[2];
rz(1.4477504) q[2];
rz(-2.0041806) q[3];
sx q[3];
rz(-1.3237938) q[3];
sx q[3];
rz(2.494273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85957134) q[0];
sx q[0];
rz(-0.30680007) q[0];
sx q[0];
rz(-0.71722537) q[0];
rz(1.2099129) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(2.4338914) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4004138) q[0];
sx q[0];
rz(-2.3614863) q[0];
sx q[0];
rz(-2.0692503) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3466481) q[2];
sx q[2];
rz(-1.7094311) q[2];
sx q[2];
rz(-0.56425205) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3860491) q[1];
sx q[1];
rz(-1.2554597) q[1];
sx q[1];
rz(2.3029033) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95146146) q[3];
sx q[3];
rz(-1.2776432) q[3];
sx q[3];
rz(0.64630634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5132961) q[2];
sx q[2];
rz(-0.6568903) q[2];
sx q[2];
rz(-0.90325242) q[2];
rz(-1.603027) q[3];
sx q[3];
rz(-0.86849803) q[3];
sx q[3];
rz(-0.85047754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83508867) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(-2.3256336) q[1];
sx q[1];
rz(-2.7201256) q[1];
sx q[1];
rz(1.0526007) q[1];
rz(-1.1076526) q[2];
sx q[2];
rz(-1.1504428) q[2];
sx q[2];
rz(3.0124315) q[2];
rz(1.611931) q[3];
sx q[3];
rz(-1.516468) q[3];
sx q[3];
rz(-0.80249912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
