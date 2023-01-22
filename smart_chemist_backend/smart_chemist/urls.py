"""smart_chemist URL definitions"""
from django.urls import path
from rest_framework.routers import DefaultRouter
from . import views

urlpatterns = [
    path('api/names', views.SmartChemistView.as_view())
]

# router = DefaultRouter()
# router.register('jobs', views.DoGSiteJobViewSet)
# router.register('info', views.DoGSiteInfoViewSet)
# urlpatterns.extend(router.urls)